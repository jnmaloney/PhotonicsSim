// Compile with 
// gcc main.cpp -lstdc++ -lm
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>


// Input
double I = 0;
double dt = 1e-9;
double t = 0;
double get_I(double t0)
{
  double h0 = 1.0;
  int blink_freq = 1.5e6;
  return h0 * (((int)(blink_freq * t0)) % 2);
}


// To calculate
double N_p = 0; // photon number
double N = 0; // carrier number


// Constants
double c = 299792458;
double q = 1.602e-19;
double h = 6.62607015e-34;

// Parameters
double eta_i = 0.75;
double lambda = 1310e-9; // Laser Wavelength
double nu = c / lambda; // Laser Frequency
double V = 4e-18; // Volume of active region
//double V_p = 0.0; // Photon Volume?
double a = 0.25e-20; // differential gain
//double a = 0.50e-20;
double epsilon = 0.7e-22;
//double epsilon = 0;
double tau = 2.0e-9;
double n_g = 4.2; // group effective index
double nu_g = c / n_g; // group velocity
double N_tr = 1.8e24; // Trancparency Number
double Gamma = 0.03; // Confinement factor

double beta_sp = 0.869e-4;
double R1 = 0.32;
double R2 = 0.32;
double L = 250.0e-6;

double alpha_i = 500.0;   // Internal loss
double alpha_m = 1.0 / (2 * L) * log(1 / (R1 * R2)); // Mirror loss
double inv_tau_sp = 1.0 / 9e-7;
double inv_tau_p = nu_g * (alpha_i + alpha_m); // Photon lifetime

//confinement factor Gamma is defined as V/V_p
double V_p = V / Gamma; // mode volume

// Find out the analytcal expected thresholds
double N_th = N_tr + inv_tau_p * (1.0 / (nu_g * a * Gamma));
double I_th = N_th * q * V / (tau * eta_i); 
double I_tr = N_tr * q * V / (tau * eta_i); 
double dPdI = eta_i * h * nu / q;
double t_d = tau * log(2.0);

// // Simulation init
// double c1 = eta_i / (q * V);
// double c2 = -1.0 / tau;
// double c3 = beta_sp * inv_tau_sp * N;
// double c4 = -beta_sp * inv_tau_p * N_p;


// double g_th = (1.0 / Gamma) * (alpha_i + alpha_m);
double step();
void init(double);


int main(int, char**)
{
  // Analytic solutions
  std::cout << "N_tr: " << N_th << std::endl;
  std::cout << "N_th: " << N_th << std::endl;
  std::cout << "I_tr: " << I_tr * 1000 << " mA" << std::endl;
  std::cout << "I_th: " << I_th * 1000 << " mA" << std::endl;
  std::cout << "slope efficiency:  " << dPdI << std::endl;
  std::cout << "t_d:  " << t_d * 1000 << " ms" << std::endl;

  // Calculate time series
  init(0);
  int k = 0;
  std::ofstream file1;
  file1.open("a.txt");
  while (t < 0.002e-3)
  {
    double P_out = step();
    ++k;
    if (k == 100)
    {
      file1 << "t: " << 1000 * t << " ms    I: " << I * 1000.0 << " mA    P_out: " << 1000.0 * P_out << " mW    N: " << N << std::endl;
      k = 0;
    }
  }
  return 0;

  // Calculate P vs I
  for (int i = 1; i < 10; ++i)
  {
    double sim_I = 0.014 + (double)i / 10.0 * 0.002;
    //double sim_I = (double)i / 100.0 * 0.036;
    //double sim_I = 0.036;
    init(sim_I);

    // Find steady state of P_out
    double eps = 1e-16;
    double P_out = 0.0;
    double dP_out = 0.0;
    do
    {
      //double old_P_out = P_out;
      double new_P_out = step();
      dP_out = fabs(P_out - new_P_out);
      P_out = new_P_out;      
      //std::cout << "P_out: " << P_out << " W" << std::endl;
      //std::cout << (fabs(dP_out < eps)) << std::endl;
      if (dP_out < eps)
      {
        //std::cout << "P_out: " << P_out << " W       dP_out: " << dP_out << " old:  " << old_P_out << std::endl;
        ++k;
      }
      else
      {
        //if (k) std::cout << "reset: " << dP_out << " W" << std::endl;
        k = 0;
      }
    } 
    while (k < 1000);
    //while (P_out == 0 || fabs(dP_out) > eps);
    //while (P_out < 1);//
    std::cout << "I:     " << sim_I * 1000 << " mA    P_out: " << 1000.0 * P_out << " mW" << std::endl;
    //std::cout << "dP_out:" << fabs(dP_out) << " W" << std::endl;
    // std::cout << (P_out == 0 || fabs(dP_out) > eps) << std::endl;
    // std::cout << (P_out == 0) << std::endl;
    // std::cout << (fabs(dP_out) > eps) << std::endl;
  }

  return 0;
}


void init(double sim_I)
{
  I = sim_I;
  t = 0;
  //dt = 1e-13;
  dt = 1e-12;
}


double c1 = eta_i / (q * V);
double c2 = -1.0 / tau;
double c3 = beta_sp * inv_tau_sp;// * N;
double c4 = -beta_sp * inv_tau_p;// * N_p;
double step()
{
  // Simulation step
  double x = nu_g * a * (N_p / (1.0 + epsilon * N_p) * (N - N_tr));
  double dN = dt * (c1 * I + c2 * N - x);
  //double dN_p = dt * Gamma * (c3 * N + c4 * N_p + x);
  double dN_p = dt * (Gamma * (nu_g * a * (N - N_tr) / (1 + epsilon * N_p) * N_p + beta_sp * N * inv_tau_sp) - N_p * inv_tau_p);

  //std::cout << c1 * I << " " << c2 * N << " " << -x << std::endl;

  // update
  t += dt;
  N += dN;
  N_p += dN_p;

  I = 0.009 * get_I(t);
  //I = 0.021 * get_I(t);
  //I = 0.0296 * get_I(t);
  //I = 0.060 * get_I(t);
  //I = 0.120 * get_I(t);

  //return dN_p;
      // std::cout << "I: " << I << "" << std::endl;
      // std::cout << "t: " << t << "" << std::endl;
      // std::cout << "dN: " << dN << "" << std::endl;
      // std::cout << "dN_p: " << dN_p << "" << std::endl;
    //std::cout << "N:   " << N << "" << std::endl;
    //std::cout << "N_p: " << N_p << "" << std::endl;

  double P_out = N_p * V_p * h * nu * alpha_m * nu_g;
  return P_out;
}


// double step2()
// {
//   // output

//   // Simulation step
//   I = get_I(t);
//   double x = nu_g * a * (N_p / (1.0 + epsilon * N_p) * (N - N_tr));
//   double dN = dt * (c1 * I + c2 * N - x);
//   double dN_p = dt * Gamma * (c3 * N + c4 * N_p + x);

//   // update
//   t += dt;
//   N += dN;
//   //N_p += dN_p;
//   return dN_p;
// }