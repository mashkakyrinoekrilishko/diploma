#include <fstream>
#include <iostream>
#include <vector>
#include <set>

using namespace std;

struct node
{
   double r;
   double z;
};

struct material
{
   double mu;
   double lambda;
};

struct element
{
   vector<int> num; //локальные узлы
   int mater;
};

struct layer {
   double x0, x1;
   double y0, y1;
   unsigned int n_mat;
};

vector<node> nodes;         // все узлы в порядке глобальной нумерации
vector<element> elems;      // все элементы в прядке глобальной нумерации
vector<material> materials; // все материалы по индексам
vector<vector<int>> KR1;    // KR1[i][j] на j-ом узле заданы i-ые краевые 1 рода
vector<layer> layers;       //слои в среде
vector<pair<int, vector<int>>>KR2_z; // граница параллельна z на j-ом узле заданы краевые 2 рода

class KURS
{
public:
   vector<double> di, al, au, b_loc, b, time, M2, M1, temp, p, z, r, x0, x1, x2, Ar, y, L, D, U;
   vector<int> ja, ia;
   vector<vector<double>> A_loc, G_loc, M_loc, M_loc_g;
   vector<double> ax, ay; //узлы областей до разбиения
   set <double> x_s, y_s; //границы КЭ
   int N, Kel, time_number, istoc, nxw, nyw, n_layers, maxIter = 10000, iter = 0;
   double eps = 1E-15, normR = 0, normB, hr, hz, rp, lambda, zs, t0, t1, t2, nu0, nu1, nu2, istoc_r, istoc_z;

   double gamma(double r); //гамма 
   double f(double r, double z, double t); //правая часть
   double f_KR1(double r, double z, int s1_id, double t); //правая часть для перых краевых
   double f_KR2(double r, double z); //правая части для вторых краевых
   double u(double r, double z, double t);
   void Input(); //считывание всех необходимых данных
   void Generate_Net(string net_file); //генерация сетки и конечных элементов
   void Read_Layers(string obj_file); //формируем слоистую среду
   void Get_G(); //получение локальной матрицы G
   void Get_M();//получение локальных матриц M
   void Get_b(double t); // получение локального b
   void Locals(int el_id, double t); // получение локальной матрицы А
   void Assemble_Locals(int el_id); // внесение локальных в глобальную СЛАУ
   void Nonlinear(); //цикл по времени
   void Get_KR1(double t); // учет первых краевых
   void Generate_Portrait(); // генерация портрета
   void multiply(vector<vector<double>>& A, vector<double>& x, vector<double>& res);
   int  Get_Num_Layer(double x0, double y0, double x1, double y1); //определяем к какому слою относится элемент
   void Get_KR2(); //учет вторых краевых

   ///////////////////для решателя/////////////////////////
   void LOS_LU(); 
   void FactLU(vector<double>& L, vector<double>& U, vector<double>& D);
   void Direct(vector<double>& L, vector<double>& D, vector<double>& y, vector<double>& b);
   void Reverse(vector<double>& U, vector<double>& x, vector<double>& y);
   double Norm(vector<double>& x);
   double mult(const vector<double>& a, const vector<double>& b);
   void Ax(vector<double>& x, vector<double>& y);
   void ATx(vector<double>& x, vector<double>& y);
};

void KURS::Generate_Net(string net_file)
{
   ifstream fin(net_file + ".txt");
   fin >> nxw;

   ax.resize(nxw);

   for (int i = 0; i < nxw; i++)
      fin >> ax[i];

   //обработка по х
   for (int j = 0; j < nxw; j++)
      x_s.insert(ax[j]);

   //read y
   fin >> nyw;
   ay.resize(nyw);

   for (int i = 0; i < nyw; i++)
      fin >> ay[i];

   //обработка по y
   for (int j = 0; j < nyw; j++)
         y_s.insert(ay[j]);

   Kel = (nxw - 1) * (nyw - 1); //n
   N = nxw * nyw; //m
   nodes.resize(N);
   // Получаем узлы
   int k = 0;
   for (set<double>::iterator j = y_s.begin(); j != y_s.end(); j++)
      for (set<double>::iterator i = x_s.begin(); i != x_s.end(); i++)
      {
         nodes[k].r = *i;
         nodes[k].z = *j;
         if (nodes[k].r == istoc_r && nodes[k].z == istoc_z)
             istoc = k;
         k++;

      }

   // Получаем КЭ
   int ku = 0;
   elems.resize(Kel);
   for (int i = 0; i < Kel; i++)
      elems[i].num.resize(4);
   set<double>::iterator j_it = y_s.begin();
   for (unsigned int j = 0; j < nyw - 1; j++)
   {
      set<double>::iterator i_it = x_s.begin();
      for (unsigned int i = 0; i < nxw - 1; i++)
      {
         elems[ku].num[0] = nxw * j + i;
         elems[ku].num[1] = nxw * j + i + 1;
         elems[ku].num[2] = nxw * (j + 1) + i;
         elems[ku].num[3] = nxw * (j + 1) + i + 1;
         double x0 = *i_it, y0 = *j_it;
         set<double>::iterator k = j_it;
         k++;
         i_it++;
         double x1 = *i_it, y1 = *k;
         unsigned int num_area = Get_Num_Layer(x0, y0, x1, y1);
         elems[ku].mater = layers[num_area].n_mat;
         ku++;
      }
      j_it++;
   }
}

int KURS::Get_Num_Layer(double x0, double y0, double x1, double y1)
{
   for (int i = 0; i < n_layers; i++)
      if (layers[i].x0 <= x0 && layers[i].x1 >= x1 && layers[i].y0 <= y0 && layers[i].y1 >= y1)
         return i;
   return 0;
}

void KURS::Read_Layers(string obj_file)
{
   ifstream fin(obj_file + ".txt");
   layer tmp;
   fin >> n_layers;
   layers.resize(n_layers);
   for (int i = 0; i < n_layers; i++) {
      fin >> tmp.x0 >> tmp.x1 >> tmp.y0 >> tmp.y1 >> tmp.n_mat;
      layers.push_back(tmp);
   }
}

double KURS::gamma(double r) // значение гамма всегда равно 1/r^2
{
 //  return 1.0;
   return 1.0 / r / r;
}

double KURS::f(double r, double z, double t) // значение f по индексу f_id 
{
    //  return (1-lambda)/r;
   //   return 10*t;
     // return -15*t+5*r*r;
      return 0;
     // return 100 * gamma(r);
      //return r * r + z * z - 6;
    //  return t * gamma(r) + 1;
}

double  KURS::f_KR1(double r, double z, int s1_id, double t) // значение краевого S1 по индексу f_id
{
   switch (s1_id)
   {
   case 0:
    //  return t + r;
      return  1./r;
     // return 100;
      //return r * r + z * z;
   default:
      cout << "can't find s1 № " << s1_id << "\n";
      break;
   }
}

double KURS::f_KR2(double r, double z) // значение краевого S1 по индексу f_id
{
        //  return t + r;
        return  0.;
        // return 100;
         //return r * r + z * z;
}

//Функция, которая является решением уравнения
double KURS::u(double r, double z, double t)
{
  // return  t + r; 
   //return 100;
   return 1./r;
   //return x*x;
}

void KURS::Input() // чтение данных
{
   ifstream in;

   in.open("istoc.txt");
   in >> istoc_r >> istoc_z;
   in.close();

   Read_Layers("layers");
   Generate_Net("net");
   int Nmat, NKR1, NKR2;
   in.open("KR1.txt");
   in >> NKR1;
   KR1.resize(NKR1);
   for (int i = 0; i < NKR1; i++)
   {
      int size;
      in >> size;
      KR1[i].resize(size);
      for (int j = 0; j < size; j++)
         in >> KR1[i][j];
   }
   in.close();

   in.open("KR2_z.txt");
   in >> NKR2;
   KR2_z.resize(NKR2);
   for (int i = 0; i < NKR2; i++)
   {
       int size;
       in >> size >> KR2_z[i].first;
       KR2_z[i].second.resize(size);
       for (int j = 0; j < size; j++)
           in >> KR2_z[i].second[j];
   }
   in.close();

   in.open("material.txt");
   in >> Nmat;
   materials.resize(Nmat);
   for (int i = 0; i < Nmat; i++)
   {
      in >> materials[i].mu;
      materials[i].lambda = 1. / materials[i].mu;
   }
   in.close();

   in.open("time.txt");
   in >> time_number; //считываем количество временных слоев
   time.resize(time_number);
   for (int i = 0; i < time_number; i++) //считываем временные слои
      in >> time[i];
}

void KURS::Get_G() // получение локальной G
{
   double a1 = (lambda * hz * rp) / (6 * hr),
      a2 = (lambda * hz) / (12),
      a3 = (lambda * hr * rp) / (6 * hz),
      a4 = (lambda * hr * hr) / (12 * hz);
   G_loc[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
   G_loc[0][1] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
   G_loc[0][2] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
   G_loc[0][3] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;

   G_loc[1][0] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
   G_loc[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
   G_loc[1][2] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
   G_loc[1][3] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;

   G_loc[2][0] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
   G_loc[2][1] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
   G_loc[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
   G_loc[2][3] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;

   G_loc[3][0] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
   G_loc[3][1] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
   G_loc[3][2] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
   G_loc[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
}

void KURS::Get_M() // прибавление локальной М
{
   double g1 = gamma(rp), g2 = gamma(rp + hr);
   M_loc_g[0][0] = M_loc_g[2][2] = hr * hz / 4 *(
      g1 * (rp / 3 + hr / 15 ) +
      g2 * (rp / 9 + 2*hr / 45 ));
   M_loc_g[0][2] = M_loc_g[2][0] = hr * hz / 12 * (
      g1 * (rp / 2 + hr / 10 ) +
      g2 * (rp / 6 + hr / 15 ));
   M_loc_g[0][3] = M_loc_g[3][0]  = hr * hz / 12 * (
      g1 * (rp / 6 + hr / 15) +
      g2 * (rp / 6 + hr / 10));
   M_loc_g[1][0] = M_loc_g[3][2] = M_loc_g[2][3] = M_loc_g[0][1] = 
      hr * hz / 4 * (
      g1 * (rp / 9 + 2 * hr / 45) +
      g2 * (rp / 9 + hr / 15));
   M_loc_g[1][1] = M_loc_g[3][3] = hr * hz / 4 * (
      g1 * (rp / 9 + hr / 15) +
      g2 * (rp / 3 + 4*hr / 15));
   M_loc_g[1][2] = M_loc_g[2][1] = hr * hz / 12 * (
      g1 * (rp / 6 + hr / 15) +
      g2 * (rp / 6 + hr / 10));
   M_loc_g[1][3] = M_loc_g[3][1] = hr * hz / 12 * (
      g1 * (rp / 6 + hr / 10) +
      g2 * (rp / 2 + 2*hr / 5));

   M_loc[0][0] = M_loc[2][2] = hr * hz / 4 * (4*rp / 9 + hr / 9);
   M_loc[0][2] = M_loc[2][0] = hr * hz / 12 * (2*rp / 3 + hr / 6);
   M_loc[0][3] = M_loc[3][0] = hr * hz / 12 * (rp / 3 + hr / 6);
   M_loc[1][0] = M_loc[3][2] = M_loc[2][3] = M_loc[0][1] = hr * hz / 4 * (2*rp / 9 + hr / 9);
   M_loc[1][1] = M_loc[3][3] = hr * hz / 4 * (4*rp / 9 + hr / 3);
   M_loc[1][2] = M_loc[2][1] = hr * hz / 12 * (rp / 3 + hr / 6);
   M_loc[1][3] = M_loc[3][1] = hr * hz / 12 * (2*rp / 3 + hr / 2);
}

void KURS::Get_b(double t) // получение локального b
{
   double f1 = f(rp, zs, t),
      f2 = f(rp + hr, zs, t),
      f3 = f(rp, zs + hz, t),
      f4 = f(rp + hr, zs + hz, t);
   b_loc[0] =
      f1 * (hr * hz / 3 * (rp / 3 + hr / 12)) +
      f2 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
      f3 * (hr * hz / 6 * (rp / 3 + hr / 12)) +
      f4 * (hr * hz / 6 * (rp / 6 + hr / 12));
   b_loc[1] =
      f1 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
      f2 * (hr * hz / 3 * (rp / 3 + hr / 4)) +
      f3 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
      f4 * (hr * hz / 6 * (rp / 3 + hr / 4));
   b_loc[2] =
      f1 * (hr * hz / 6 * (rp / 3 + hr / 12)) +
      f2 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
      f3 * (hr * hz / 3 * (rp / 3 + hr / 12)) +
      f4 * (hr * hz / 3 * (rp / 6 + hr / 12));
   b_loc[3] =
      f1 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
      f2 * (hr * hz / 6 * (rp / 3 + hr / 4)) +
      f3 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
      f4 * (hr * hz / 3 * (rp / 3 + hr / 4));
}

void KURS::Locals(int el_id, double t) // получение локальной матрицы А
{
   element el = elems[el_id];
   hr = nodes[el.num[1]].r - nodes[el.num[0]].r;
   hz = nodes[el.num[2]].z - nodes[el.num[0]].z;
   rp = nodes[el.num[0]].r;
   zs = nodes[el.num[0]].z;
   lambda = materials[el.mater].lambda;
   Get_G(); 
   Get_M();
   Get_b(t);
   //for (int i = 0; i < 4; i++)
   //   for (int j = 0; j < 4; j++)
   //    //  A_loc[i][j] = G_loc[i][j] + M_loc[i][j] * nu0;
   //      A_loc[i][j] = G_loc[i][j] + M_loc[i][j]*nu0 + M_loc_g[i][j];
}

void KURS::Generate_Portrait() // генерация портрета
{
   ia.resize(N + 1);
   ja.resize(16 * Kel);
   vector<int> t1(16 * Kel), t2(16 * Kel), beg(N);
   int s = 0;
   for (int i = 0; i < N; i++)
      beg[i] = 0;
   for (int el = 0; el < Kel; el++)
   {
      for (int i = 0; i < 4; i++) 
      {
         int k = elems[el].num[i];
         for (int j = i + 1; j < 4; j++)
         {
            int ind1 = k;
            int ind2 = elems[el].num[j];
            if (ind2 < ind1)
            {
               ind1 = ind2;
               ind2 = k;
            }
            int iaddr = beg[ind2];
            if (iaddr == 0)
            {
               s++;
               beg[ind2] = s;
               t1[s] = ind1;
               t2[s] = 0;
            }
            else
            {
               while (t1[iaddr] < ind1 && t2[iaddr] > 0)
                  iaddr = t2[iaddr];
               if (t1[iaddr] > ind1)
               {
                  s++;
                  t1[s] = t1[iaddr];
                  t2[s] = t2[iaddr];
                  t1[iaddr] = ind1;
                  t2[iaddr] = s;
               }
               else if (t1[iaddr] < ind1)
               {
                  s++;
                  t2[iaddr] = s;
                  t1[s] = ind1;
                  t2[s] = 0;
               }
            }
         }
      }
   }

   ia[0] = 0;
   for (int i = 0; i < N; i++)
   {
      ia[i + 1] = ia[i];
      int a = beg[i];
      while (a != 0)
      {
         ja[ia[i + 1]] = t1[a];
         ia[i + 1]++;
         a = t2[a];
      }
   }
}

void KURS::Assemble_Locals(int el_id) // внесение локальных A, b  в глобальную СЛАУ
{
   vector<int> L = elems[el_id].num;
   int k = elems[el_id].num.size(); // размерность локальной матрицы
   for (int i = 0; i < k; i++)
      di[L[i]] += A_loc[i][i];

   for (int i = 0; i < 4; i++)
      for (int j = 0; j < i; j++)
         for (int k = ia[L[i]]; k < ia[L[i] + 1]; k++)
            if (ja[k] == L[j])
            {
               al[k] += A_loc[i][j];
               au[k] += A_loc[j][i];
               k++;
               break;
            }
}

// Умножение матрицы A на вектор x
// Результат в res
void KURS :: multiply(vector<vector<double>>& A, vector<double>& x, vector<double>& res)
{
   for (int i = 0; i < 4; i++)
   {
      res[i] = 0.;
      for (int j = 0; j < 4; j++)
         res[i] += A[i][j] * x[j];
   }
}

void KURS::Nonlinear()
{
   x0.resize(N);
   x1.resize(N);
   x2.resize(N);
   temp.resize(N);

   //Инициализация решения на 0-м слое как решение стационарной задачи
   au.clear();
   au.resize(ia[N]);
   al.clear();
   al.resize(ia[N]);
   di.clear();
   di.resize(N);
   b.clear();
   b.resize(N);
   ja.resize(ia[N]);
   G_loc.resize(4);
   M_loc.resize(4);
   A_loc.resize(4);
   b_loc.resize(4);
   M_loc_g.resize(4);
   for (int i = 0; i < 4; i++)
   {
       A_loc[i].resize(4);
       M_loc[i].resize(4);
       M_loc_g[i].resize(4);
       G_loc[i].resize(4);
   }
   for (int i = 0; i < Kel; i++)
   {
       Locals(i, time[0]);
       for (int i = 0; i < 4; i++)
           for (int j = 0; j < 4; j++)
               A_loc[i][j] = G_loc[i][j] + M_loc_g[i][j];
       Assemble_Locals(i); //сборка левой части
       //сборка вектора правой части
       for (int j = 0; j < 4; j++)
           b[elems[i].num[j]] += b_loc[j];
   }
   b[istoc] = b[istoc] + 10e10 / 2 / 3.14 / nodes[istoc].r; //задаем точечный источник в узле сетки
   Get_KR2();
   Get_KR1(time[0]);
   LOS_LU();
   ofstream out("result.txt");
   out << "t = " << time[0] << endl;
   out << "Полученное решение:" << endl;
   for (int j = N - 1; j >= 0; j--)
   {
       out << x0[j] << "\t";
       if (j % 21 == 0 && j != (N - 1))
           out << endl;
   }
   out << endl;
   for (int i = 0; i < N; i++) //запоминаем новые старые вектора
   {
       x2[i] = x0[i];
       x0[i] = 0.;
   }
   //////////////////////////////////////////////////////////////////

   //Инициализация решения на 1-м слое
   au.clear();
   au.resize(ia[N]);
   al.clear();
   al.resize(ia[N]);
   di.clear();
   di.resize(N);
   b.clear();
   b.resize(N);
   ja.resize(ia[N]);
   G_loc.resize(4);
   M_loc.resize(4);
   A_loc.resize(4);
   b_loc.resize(4);
   M_loc_g.resize(4);
   for (int i = 0; i < 4; i++)
   {
       A_loc[i].resize(4);
       M_loc[i].resize(4);
       M_loc_g[i].resize(4);
       G_loc[i].resize(4);
   }
   M2.resize(N);
   //считаем коэффициенты по времени для двухслойной неявной схемы
   t0 = time[1] - time[0]; //t0
   nu0 = 1 / t0;
   for (int i = 0; i < Kel; i++)
   {
       Locals(i, time[1]);
       for (int i = 0; i < 4; i++)
           for (int j = 0; j < 4; j++)
               A_loc[i][j] = G_loc[i][j] +M_loc[i][j] * nu0 + M_loc_g[i][j];
       Assemble_Locals(i); //сборка левой части
       //сборка вектора правой части
       for (int j = 0; j < 4; j++)
           temp[j] = x2[elems[i].num[j]];
       multiply(M_loc, temp, M2);//M*x2
       for (int j = 0; j < 4; j++)
           M2[j] = M2[j] * nu0;
       for (int j = 0; j < 4; j++)
           b[elems[i].num[j]] += b_loc[j] + M2[j];
   }
//   b[istoc] = b[istoc] + 10e10 / 2 / 3.14 / nodes[istoc].r; //не задаем точечный источник в узле сетки
   Get_KR2();
   Get_KR1(time[1]);
   LOS_LU();
   out << "t = " << time[1] << endl;
   out << "Полученное решение:" << endl;
   for (int j = N - 1; j >= 0; j--)
   {
       out << x0[j] << "\t";
       if (j % 21 == 0 && j != (N - 1))
           out << endl;
   }
   out << endl;
   for (int i = 0; i < N; i++) //запоминаем новые старые вектора
   {
       x1[i] = x0[i];
       x0[i] = 0.;
   }
   //////////////////////////////////////////////////////////////////////////////////////////////

  /* for (int i = 0; i < N; i++)
   {
      if (i == istoc)
      {
         x2[i] = -10e10 / 2 / 3.14 / nodes[istoc].r / 3 * nodes[i].r * (2 - nodes[i].r);
         x1[i] = -10e10 / 2 / 3.14 / nodes[istoc].r / 3 * nodes[i].r * (2 - nodes[i].r);
      }
      else
      {
         x2[i] = u(nodes[i].r, nodes[i].z, time[0]);
         x1[i] = u(nodes[i].r, nodes[i].z, time[1]);
      }
   }*/
   M1.resize(N);
   M2.resize(N);
 //  ofstream out("result.txt");
   for (int t = 2; t < time_number; t++)	//цикл по временным слоям, начиная с 3-го
   {
      //очищаем вектора для каждого временного слоя
      au.clear();
      au.resize(ia[N]);
      al.clear();
      al.resize(ia[N]);
      di.clear();
      di.resize(N);
      b.clear();
      b.resize(N);
      ja.resize(ia[N]);
      G_loc.resize(4);
      M_loc.resize(4);
      A_loc.resize(4);
      b_loc.resize(4);
      M_loc_g.resize(4);
      for (int i = 0; i < 4; i++)
      {
         A_loc[i].resize(4);
         M_loc[i].resize(4);
         M_loc_g[i].resize(4);
         G_loc[i].resize(4);
      }
      //считаем коэффициенты по времени
      t0 = time[t] - time[t - 1]; //t0
      t1 = time[t] - time[t - 2]; //t
      t2 = time[t - 1] - time[t - 2]; //t1
      nu0 = (t1+t0)/t1/t0;
      nu1 = t1/t2/t0;
      nu2 = t0/t2/t1;
      for (int i = 0; i < Kel; i++)
      {
         Locals(i, time[t]);
         for (int i = 0; i < 4; i++)
             for (int j = 0; j < 4; j++)
                 A_loc[i][j] = G_loc[i][j] + M_loc[i][j] * nu0 + M_loc_g[i][j];
         Assemble_Locals(i); //сборка левой части
         //сборка вектора правой части
         //F = b_loc(j)-nu2*M*x2+nu1*M*x1
         for (int j = 0; j < 4; j++)
            temp[j] = x2[elems[i].num[j]];
         multiply(M_loc, temp, M2);//M*x2
         for (int j = 0; j < 4; j++)
            temp[j] = x1[elems[i].num[j]];
         multiply(M_loc, temp, M1);//M*x1
         for (int j = 0; j < 4; j++)
         {
            M1[j] = M1[j] * nu1;
            M2[j] = M2[j] * nu2;
         }
         for (int j = 0; j < 4; j++)
            b[elems[i].num[j]] += b_loc[j] - M2[j] + M1[j];
      }
  //    b[istoc] = b[istoc] + 10e10/2/3.14/nodes[istoc].r; //задаем точечный источник в узле сетки
      Get_KR2();
      Get_KR1(time[t]);
      LOS_LU();
      out << "t = " << time[t] << endl;
      //out << "Аналитическое решение:" << endl;
      //for (int j = 0; j < N; j++)
      //   if (j == istoc)
      //      temp[j] = -10e6 / 2 / 3.14 / nodes[istoc].r/3*nodes[j].r*(2- nodes[j].r);
      //   else 
      //    temp[j] = u(nodes[j].r, nodes[j].z, time[t]);
      //for (int j = N-1; j >= 0; j--)
      //{
      //   out << temp[j] << "\t";
      //   if (j % 10 == 0 && j!=(N-1))
      //      out << endl;
      //}
  //    out << "Полученное решение:" << endl;
      for (int j = N - 1; j >= 0; j--)
      {
         out << x0[j] << "\t";
         if (j % 21 == 0 && j != (N - 1))
            out << endl;
      }
     /* out << "Абсолютная погрешность:" << endl;
      for (int j = N - 1; j >= 0; j--)
      {
         out << temp[j] - x0[j] << "\t";
         if (j % 10 == 0 && j != (N - 1))
            out << endl;
      }*/
      out << endl ;
      for (int i = 0; i < N; i++) //запоминаем новые старые вектора
      {
         x2[i] = x1[i];
         x1[i] = x0[i];
         x0[i] = 0.;
      }
   }
}

void KURS::Get_KR1(double t) // учет первых краевых
{
   for (int i = 0; i < KR1.size(); i++)
      for (int j = 0; j < KR1[i].size(); j++)
      {
         int node_id = KR1[i][j];
         di[node_id] = 1;
         b[node_id] = f_KR1(nodes[node_id].r, nodes[node_id].z, i, t);
         for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
            al[k] = 0;
         for (int k = 0; k < ja.size(); k++)
            if (ja[k] == node_id)
               au[k] = 0;
      }
}

void KURS::Get_KR2() // учет вторых краевых MS.b - вектор вклада от краевых
{
    for (int i = 0; i < KR2_z.size(); i++)
        for (int j = 0; j < KR2_z[i].second.size() - 1; j++)
        {
            int node_id1 = KR2_z[i].second[j],node_id2 = KR2_z[i].second[j + 1];
            double hz = nodes[node_id2].z - nodes[node_id1].z,
                f1 = f_KR2(nodes[node_id1].r, nodes[node_id1].z),
                f2 = f_KR2(nodes[node_id1].r, nodes[node_id1].z + hz);
            b_loc[0] = f1 * (hz / 3) + f2 * (hz / 6);
            b_loc[1] = f1 * (hz / 6) + f2 * (hz / 3);
            b[node_id1] += b_loc[0];
            b[node_id2] += b_loc[1];
        }
}

////////////////решатель///////////////////////////////////
void KURS::ATx(vector<double>& x, vector<double>& y)
{
   for (int i = 0; i < N; i++)
   {
      y[i] = di[i] * x[i];
      for (int j = ia[i]; j < ia[i + 1]; j++)
      {
         int k = ja[j];
         y[i] += au[j] * x[k];
         y[k] += al[j] * x[i];
      }
   }
}

void KURS::Ax(vector<double>& x, vector<double>& y)
{
   for (int i = 0; i < N; i++)
   {
      y[i] = di[i] * x[i];
      for (int j = ia[i]; j < ia[i + 1]; j++)
      {
         int k = ja[j];
         y[i] += al[j] * x[k];
         y[k] += au[j] * x[i];
      }
   }
}

double KURS::Norm(vector<double>& x)
{
   double norm = 0;
   for (int i = 0; i < N; i++)
      norm += x[i] * x[i];
   return(sqrt(norm));
}

double KURS::mult (const vector<double>& a, const vector<double>& b)
{
   double res = 0;
   for (int i = 0; i < a.size(); i++)
      res += a[i] * b[i];
   return res;
}

void KURS::LOS_LU()
{
   r.resize(N);
   z.resize(N);
   p.resize(N);
   Ar.resize(N);
   y.resize(N);
   L.resize(ia[N]);
   D.resize(N);
   U.resize(ia[N]);
   normB = Norm(b);
   cout.precision(15);
   FactLU(L, U, D);
   double p_p = 0, p_r = 0, r_r = 0, Ar_p = 0;
   double a = 0, B = 0, eps2 = 1e-10;
   Ax(x0, y); // y = A * x0
   // y = B - A * x0
   for (int i = 0; i < N; i++)
      y[i] = b[i] - y[i];
   Direct(L, D, r, y); // r0 = L^(-1) * (B - A *x0)
   Reverse(U, z, r); // z0 = U^(-1) * r0
   Ax(z, y); // y = A * z0
   Direct(L, D, p, y); // p0 = L^(-1) * (A *  z0)
   r_r = mult(r, r);
   normR = sqrt(r_r) / normB;
   for (iter = 1; iter < maxIter + 1 && normR >= eps; iter++)
   {
      p_p = mult(p,p);
      p_r = mult(p,r);
      a = p_r / p_p;
      // x(k) = x(k-1) + a(k) * z(k-1)
      // r(k) = r(k-1) - a(k) * p(k-1)
      for (int i = 0; i < N; i++)
      {
         x0[i] = x0[i] + z[i] * a;
         r[i] = r[i] - p[i] * a;
      }
      Reverse(U, y, r); // y = U^(-1) *  r(k)
      Ax(y, Ar); // Ar = A * U^(-1)  *r(k)
      Direct(L, D, Ar, Ar); // Ar = L^(-1) * A *U ^ (-1)* r(k)
      Ar_p = mult(Ar,p); // (Ar, p)
      B = -(Ar_p / p_p);
      // z(k) = U^(-1) * r(k) + B(k) * z(k-1) 
      // p(k) = L^(-1) * A * U^(-1) * r(k) + B(k)* p(k - 1)
      for (int i = 0; i < N; i++)
      {
         z[i] = y[i] + z[i] * B;
         p[i] = Ar[i] + p[i] * B;
      }
      if (r_r - (r_r - a * a * p_p) < eps2)
         r_r = mult(r,r);
      else
         r_r = r_r - a * a * p_p;
      normR = sqrt(r_r) / normB;
      cout << iter << ". " << normR << endl;
   }
}

void KURS::FactLU(vector<double>& L, vector<double>& U, vector<double>& D)
{
   L = al;
   U = au;
   D = di;
   double l, u, d;
   for (int k = 0; k < N; k++)
   {
      d = 0;
      int i0 = ia[k], i1 = ia[k + 1];
      int i = i0;
      for (; i0 < i1; i0++)
      {
         l = 0;
         u = 0;
         int j0 = i, j1 = i0;
         for (; j0 < j1; j0++)
         {
            int t0 = ia[ja[i0]], t1 = ia[ja[i0] + 1];
            for (; t0 < t1; t0++)
            {
               if (ja[j0] == ja[t0])
               {
                  l += L[j0] * U[t0];
                  u += L[t0] * U[j0];
               }
            }
         }
         L[i0] -= l;
         U[i0] -= u;
         U[i0] /= D[ja[i0]];
         d += L[i0] * U[i0];
      }
      D[k] -= d;
   }
}

// L*y = B
void KURS::Direct(vector<double>& L, vector<double>& D, vector<double>& y, vector<double>& b)
{
   y = b;
   for (int i = 0; i < N; i++)
   {
      double sum = 0;
      int k0 = ia[i], k1 = ia[i + 1];
      int j;
      for (; k0 < k1; k0++)
      {
         j = ja[k0];
         sum += y[j] * L[k0];
      }
      double buf = y[i] - sum;
      y[i] = buf / D[i];
   }
}

// U*x = y
void KURS::Reverse(vector<double>& U, vector<double>& x, vector<double>& y)
{
   x = y;
   for (int i = N - 1; i >= 0; i--)
   {
      int k0 = ia[i], k1 = ia[i + 1];
      int j;
      for (; k0 < k1; k0++)
      {
         j = ja[k0];
         x[j] -= x[i] * U[k0];
      }
   }
}

int main()
{
   KURS A;
   A.Input();
   A.Generate_Portrait();
   A.Nonlinear();
   return 0;
}
