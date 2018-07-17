
#define _USE_MATH_DEFINES
#include <math.h>

#include "Units.h"
#include <typeinfo>
#include <iostream>
#include <cassert>

using namespace Units;

void mechanic()
{
   Length l;
   Velocity v;
   Time t;
   Energy E;
   Acceleration a;
   Mass m;
   Force F;
   Frequency f;
   Stiffness k;
   MassDensity r;

   m = r * l * l *l;

   // free motion
   l = v * t;
   
   //Accelerated motion
   v = a * t;
   l = a * t * t / 2;
   t = sqrt(2 * l / a);
   
   //Kinetic energy
   E = m * v * v / 2;
   v = sqrt(E / m);

   //Newton law
   F = m * a;

   // Gravitational acceleration
   F = m * g;

   // work of force
   E = F * l;

   // Newton's law of universal gravitation
   F = (m * m) / (l * l) * G;

   // Spring pendulum frequency
   f = sqrt(k / m);

   //string frequency
   f = sqrt(F / r) / l / l;
}

void electrostatic()
{
   Charge Q = 1.0_C;
   Length l;
   Force F;
   ElStrength E;
   ElDipoleMoment p;
   ElectricFlux Fi;

   //Gauss law
   Fi = Q / e0;
   E = Q / l / l / e0;

   //Coulomb's law
   F = Q * Q / ( l * l * e0);

   //Dipole field
   E = p / (e0 * l * l * l);
}

void termodynamics()
{
   Energy E;
   Temperature T;
   Mass m;
   Pressure p;
   Density n;
   Volume V;
   SubstanceAmount N;
   Velocity v;
   MolarMass nu;
   Force F;
   DynViscosity S;
   KinViscosity vis;
   Length l;
   Area s;
   MassDensity ro;
   Charge Q;

   E = T * kb;
   E = kb * T;
   T = E / kb;
   SpecificHeat c;
   E = c * m * T;
 
   //
   p = kb * T * n;

   //Boyle's law
   Pressure p2;
   Volume V2;
   V = p2 / p * V2;

   //Gay-Lussac's law
   Temperature T2;
   p = p2 / T2 * T;

   // Ideal gas law
   p * V ==  R * T * N;

   // sound velocity in gas
   v = sqrt(R * T / nu);

   //Dynamic viscosity
   F = S * s * v / l;
   vis = S / ro;

   //Faraday's laws of electrolysis
   m = (Q / Fc) * nu;
}

void quantumMechanics()
{
   Energy E;
   Frequency f;
   Length l;
   Impulse p;
   WaveNumber k;
   Time t;

   // Quantum energy
   E = h * f;
   f = E / h;

   // Matter wave
   l = h / p;
   p = h / l;
   p = h * k;
   k = p / h;

   // Uncertainty principle
   l * p <= h;
   t * E <= h;

}

void elecricity()
{
   ElPotential U = 1.0_V;
   ElCurrent I = 1.0_A;
   ElResistance R = 1.0_Ohm;
   ElConductance S = 1.0_S;
   Energy E = 1.0_J;
   Power W = 1.0_W;
   Charge Q = 1.0_C;
   ElCapacity C = 1.0_F;
   ElInductance L;
   Frequency f;
   Time t;
   Length l;
   EnergyDensity w;
   Force F;

   // Ohm law
   U = I * R;
   R = U / I;
   I = U / R;
   S = I / U;
   S = 1 / R;
   R = 1 / S;
   assert(5.0_A * 2.0_Ohm == 10.0_V);
   assert(15.0_V / 5.0_Ohm == 3.0_A);
   assert(15.0_V / 3.0_A == 5.0_Ohm);
   assert(6.0_A / 3.0_V == 2.0_S);
   assert(1 / 0.5_Ohm == 2.0_S);
   assert(1 / 0.5_S == 2.0_Ohm);

   // Joule–Lenz law
   W = U * I;
   U = sqrt(W * R);
   I = sqrt(W / R);
   assert(5.0_V * 4.0_A == 20.0_W);
   assert(sqrt(25.0_W / 1.0_Ohm) == 5.0_A);
   assert(sqrt(5.0_W * 5.0_Ohm) == 5.0_V);

   //Charge
   Q = I * t;
   assert(2.0_A * 6.0_s == 12.0_C);
   
   //Capasitor charge
   Q = C * U;
   assert(2.0_V * 9.0_F == 18.0_C);

   //Capasitor resistance
   S = C * f;
   R = t / C;

   //Charge energy
   E = Q * U;

   //Capasitor energy
   E = C * U * U / 2;

   //Capasitance
   C = e0 * l;

   //Inductivity resistance
   R = L / t;
   S = t / L;


   // Oscillatory circuit frequency
   f = sqrt(1 / C / L);
  // assert(sqrt(1 / (240.0_nF * 360.0_nH)) == 542.0_kHz);


   //Inductivity energy
   E = I * I * L / 2;

   U = Q / e0 / l;

   // Field energy density
   w = (U / l) * (U / l) * e0;
   w = U * U / ( l * l) * e0;

}

void plasmaPhysics()
{
   Charge c;
   Length l;
   Temperature T;
   Density n;
   Frequency f;
   Mass m;
   MgStrength B;

   //Debye length
   l = sqrt(kb * T * e0 / (c * c  * n));
   
   //Plasma frequency
   f = sqrt(4 * M_PI * n * c * c / m / e0 * m);

   // Cyclotron resonance
   f = c * B / m;
}

void optics()
{
   Temperature T;
   Frequency f;
   SpectralRadiance B;

   //Stefan–Boltzmann law
   auto R = T * T * T * T * s;
   
   //Planck's law
   B = 2 * h * f * f * f / (c * c) / (exp(h * f / kb / T) + 1);

   f + f;
}

void magnetic()
{
   MgStrength B;
   ElCurrent I;
   Length r;
   Force F;
   Velocity v;
   Charge q;
   Time t;
   ElPotential U;

   // Ampère's force law
   F = I * r * B;

   //Lorentz force
   F = q * v * B;

   //Coil field
   B = 0.5 * mu0 * I / r;

   //Faraday's law of induction
   U = M_PI * r * r * B / t;
}

void commonOps()
{
   Length l = 2.0_m;
   
   assert(std::string(l.getName()) == "Length");
   assert(std::string((l*l).getName()) == "Area");
   assert(std::string((l*l*l).getName()) == "Volume");
   assert(std::string((l*l*l*l).getName()) == "unnamed");

   assert(std::string(l.getUnit()) == "m");
   assert(std::string(Charge::getUnit()) == "C");

   // Comparisons
   assert(l == 2.0_m);
   assert(l != 3.0_m);
   assert(l < 3.0_m);
   assert(l > 1.0_m);
   assert(l <= 3.0_m);
   assert(l >= 1.0_m);
   assert(l <= 2.0_m);
   assert(l >= 2.0_m);

   l = l + l;
   assert(l == 4.0_m);
   l += l;
   assert(l == 8.0_m);

   l = l / 2;
   assert(l == 4.0_m);
   l /= 2;
   assert(l == 2.0_m);

   l = l * 2;
   assert(l == 4.0_m);
   l *= 2;
   assert(l == 8.0_m);
   l = 2 * l;
   assert(l == 16.0_m);
   l = l - 1.0_m;
   assert(l == 15.0_m);
   l -= 1.0_m;
   assert(l == 14.0_m);
   l = -l;
   assert(l == -14.0_m);
   l = 0.0_m - l;
   assert(l == 14.0_m);

   auto r =l / l;
}

int main()
{
   commonOps();
   elecricity();
   return 0;
}