#pragma once

#ifndef UNITS_H_
#define UNITS_H_

#include <string>
#include <math.h>

#define IN1 const Unit<l, m, t, i, k, n, j, a>& rhs
#define IN2 const Unit<L, M, T, I, K, N, J, A>& rhs
#define FIN const Value & rhs
#define TMPLT template<int l, int m, int t, int i, int k, int n, int j, int a>
#define OTHR_TMPLT template<int L, int M, int T, int I, int K, int N, int J, int A>
#define RET return typename UnitName <l, m, t, i, k, n, j, a>::Name
#define B_RET return bool
#define MLT_RET return typename UnitName <l + L, m + M, t + T, i + I, k + K, n + N, j + J, a + A>::Name
#define DIV_RET return typename UnitName <l - L, m - M, t - T, i - I, k - K, n - N, j - J, a - A>::Name
#define INV_RET return typename UnitName <-l, -m, -t, -i, -k, -n, -j, -a>::Name
#define SQRT_RET typename UnitName <l / 2, m / 2, t / 2, i / 2, k / 2, n / 2, j / 2, a / 2>::Name

#define OPERATOR(op, ret, sp, ...) \
inline auto operator op (__VA_ARGS__) sp {ret;}

#define VALUE(name, unit, ...) \
typedef Unit<__VA_ARGS__> name; \
template <> \
struct UnitName<__VA_ARGS__> \
{  \
   typedef name Name; \
   static const char* getUnit() { return #unit; } \
   static const char* getName() {return #name;} \
};

#define PR_FEMTO f
#define PR_PICO p
#define PR_NANO n
#define PR_MICRO mc
#define PR_MILLI m
#define PR_KILO k
#define PR_MEGA M
#define PR_GIGA G
#define PR_TERA T

#define LITERAL(Type, ltr, coeff)          \
inline auto operator ""_##ltr(Value value) \
{                                          \
   return Type::Name(value * coeff);       \
}

#define LITERALS(Type, ltr, coeff)         \
LITERAL(Type, f##ltr, 1e-15 * coeff)       \
LITERAL(Type, p##ltr, 1e-12 * coeff)       \
LITERAL(Type, n##ltr, 1e-9 * coeff)        \
LITERAL(Type, mc##ltr, 1e-6 * coeff)       \
LITERAL(Type, m##ltr, 1e-3 * coeff)        \
LITERAL(Type, ltr, 1 * coeff)              \
LITERAL(Type, k##ltr, 1e+3 * coeff)        \
LITERAL(Type, M##ltr, 1e+6 * coeff)        \
LITERAL(Type, G##ltr, 1e+9 * coeff)        \
LITERAL(Type, T##ltr, 1e+12* coeff)

namespace Units
{
   typedef long double Value;

   TMPLT class Unit;

   TMPLT struct UnitName
   {
      typedef Unit<l, m, t, i, k, n, j, a> Name;

      // General form of unit - dimensional formula
      // It can be spicified for some derrived SI units by VALUE macro
      static const char* getUnit()
      {
         static std::string buf =
            (l != 0) ? std::string("m^") + std::to_string(l) : "" +
            (m != 0) ? std::string("kg^") + std::to_string(m) : "" +
            (t != 0) ? std::string("s^") + std::to_string(t) : "" +
            (i != 0) ? std::string("A^") + std::to_string(i) : "" +
            (k != 0) ? std::string("K^") + std::to_string(k) : "" +
            (n != 0) ? std::string("mol^") + std::to_string(n) : "" +
            (j != 0) ? std::string("cd^") + std::to_string(j) : "";
         return buf.c_str();
      }
      static const char* getName() { return "unnamed"; }
   };
   //Dimensionless numbers
   template <>
   struct UnitName<0, 0, 0, 0, 0, 0, 0, 0> { typedef Value Name; };
   
   //Base SI Units
   VALUE(Length, m, 1, 0, 0, 0, 0, 0, 0, 0);
   VALUE(Mass, kg, 0, 1, 0, 0, 0, 0, 0, 0);
   VALUE(Time, s, 0, 0, 1, 0, 0, 0, 0, 0);
   VALUE(ElCurrent, A, 0, 0, 0, 1, 0, 0, 0, 0);
   VALUE(Temperature, K, 0, 0, 0, 0, 1, 0, 0, 0);
   VALUE(SubstanceAmount, mol, 0, 0, 0, 0, 0, 1, 0, 0);
   VALUE(LumIntensity, cd, 0, 0, 0, 0, 0, 0, 1, 0);
   VALUE(Angle, cd, 0, 0, 0, 0, 0, 0, 0, 1);
   VALUE(SolidAngle, cd, 0, 0, 0, 0, 0, 0, 0, 2);
   //SI derived units
   VALUE(Area, m^2, 2, 0, 0, 0, 0, 0, 0, 0);
   VALUE(Volume, m^3, 3, 0, 0, 0, 0, 0, 0, 0);
   VALUE(Velocity, m/s, 1, 0, -1, 0, 0, 0, 0, 0);
   VALUE(ArealVelocity, m^2/s, 2, 0, -1, 0, 0, 0, 0, 0);
   VALUE(Acceleration, m/s^2, 1, 0, -2, 0, 0, 0, 0, 0);
   VALUE(Force, N, 1, 1, -2, 0, 0, 0, 0, 0);
   VALUE(Charge, C, 0, 0, 1, 1, 0, 0, 0, 0);
   VALUE(Density, m^-3, -3, 0, 0, 0, 0, 0, 0, 0);
   VALUE(Frequency, Hz, 0, 0, -1, 0, 0, 0, 0, 0);
   VALUE(Impulse, kg*m/s, 1, 1, -1, 0, 0, 0, 0, 0);
   VALUE(Pressure, Pa, -1, 1, -2, 0, 0, 0, 0, 0);
   VALUE(Energy, J, 2, 1, -2, 0, 0, 0, 0, 0);
   VALUE(GrPotential, m^2/s^2, 2, 0, -2, 0, 0, 0, 0, 0);
   VALUE(Power, W, 2, 1, -3, 0, 0, 0, 0, 0);
   VALUE(ElPotential, V, 2, 1, -3, -1, 0, 0, 0, 0);
   VALUE(ElResistance, Omh, 2, 1, -3, -2, 0, 0, 0, 0);
   VALUE(ElResistivity, Omh*m, 3, 1, -3, -1, 0, 0, 0, 0);
   VALUE(ElConductance, S, -2, -1, 3, 2, 0, 0, 0, 0);
   VALUE(ElConductivity, S/m, -3, -1, 3, 1, 0, 0, 0, 0);
   VALUE(ElCapacity, F, -2, -1, 4, 2, 0, 0, 0, 0);
   VALUE(ElInductance, H, 2, 1, -2, -2, 0, 0, 0, 0);
   VALUE(ElDipoleMoment, C*m^2, 1, 0, 1, 1, 0, 0, 0, 0);
   VALUE(ElStrength, V/m, 1, 1, -3, -1, 0, 0, 0, 0);
   VALUE(SpecificHeat, J/(kg*k), 2, 0, -2, 0, -1, 0, 0, 0);
   VALUE(WaveNumber, 1/m, -1, 0, 0, 0, 0, 0, 0, 0);
   VALUE(Stiffness, N/m, 0, 1, -2, 0, 0, 0, 0, 0);
   VALUE(InertiaMoment, kg*m^2, 2, 1, 0, 0, 0, 0, 0, 0);
   VALUE(ChargeMassRatio, C/m, 0, -1, 0, 1, 1, 0, 0, 0);
   VALUE(MassDensity, kg/m^3, -3, 1, 0, 0, 0, 0, 0, 0);
   VALUE(ChargeDensity, C/m^3, -3, 0, 1, 1, 0, 0, 0, 0);
   VALUE(CurrentDensity, A/m^2, -2, 0, 0, 1, 0, 0, 0, 0);
   VALUE(VolumetricFlow, m^3/s, 3, 0, -1, 0, 0, 0, 0, 0);
   VALUE(DensityOfFlow, 1/(m^2*s), -2, 0, -1, 0, 0, 0, 0, 0);
   VALUE(Luminance, lx, -2, 0, 0, 0, 0, 0, 1, 0);
   VALUE(MagneticMoment, J/T, 2, 0, 0, -1, 0, 0, 0, 0);
   VALUE(StefansConst, W/(m^2*K^4), 0, 1, -3, 0, -4, 0, 0, 0);
   VALUE(AvogadroConst, 1/mol, 0, 0, 0, 0, 0, -1, 0, 0);
   VALUE(MgStrength, T, 0, 1, -2, -1, 0, 0, 0, 0);
   VALUE(RadExitance, W/m^2, 0, 1, -3, 0, 0, 0, 0, 0);
   VALUE(MolarMass, kg/mol, 0, 1, 0, 0, 0, -1, 0, 0);
   VALUE(DynViscosity, Pl, -1, 1, -1, 0, 0, 0, 0, 0);

   VALUE(FaradayConstant, s*A/mol, 0, 0, 1, 1, 0, -1, 0, 0);
   VALUE(BoltsmanConstant, J/K, 2, 1, -2, 0, -1, 0, 0, 0);
   VALUE(MgConst, V*s/(A*m), 1, 1, -2, -2, 0, 0, 0, 0);
   VALUE(PlanckConstant, J*s, 2, 1, -1, 0, 0, 0, 0, 0);
   VALUE(GrConstant, m^3/(kg*s^2), 3, -1, -2, 0, 0, 0, 0, 0);
   VALUE(GasConstant, J/(K*mol), 2, 1, -2, 0, -1, -1, 0, 0);
   VALUE(ElConst, N*m^2/C^2, -3, -1, 4, 2, 0, 0, 0, 0);

   typedef Pressure                 EnergyDensity;
   typedef Energy                   Work;
   typedef WaveNumber               LensPower;
   typedef Length                   FocalLength;
   typedef Pressure                 Stress;
   typedef WaveNumber               LinearDencity;
   typedef BoltsmanConstant         HeatCapacity;
   typedef ElResistivity            ElectricFlux;
   typedef Energy                   Torque;
   typedef Stiffness                SurfaceTension;
   typedef RadExitance              Radiosity;
   typedef Frequency                RadActivity;
   typedef Frequency                Activity;
   typedef ArealVelocity            KinViscosity;
   typedef Stiffness                SpectralRadiance;
   typedef SpecificHeat             RadAbsDose;
   typedef Pressure                 YoungMdl;

   TMPLT class Unit : public UnitName<l, m, t, i, k, n, j, a>
   {
      Value val_ = 0;
   public:
      Unit() = default;
      Unit(Value val) : val_(val) {}
      inline Value getSI() const { return val_; }

      OPERATOR(-, RET(-val_), const,)
      OTHR_TMPLT OPERATOR(*, MLT_RET(val_ * rhs.getSI()), const, IN2)
      OTHR_TMPLT OPERATOR(/, DIV_RET(val_ / rhs.getSI()), const, IN2)
      OPERATOR(*=, val_ *= rhs; return *this; ,, FIN)
      OPERATOR(/=, val_ /= rhs; return *this; ,, FIN)
      OPERATOR(+=, val_ += rhs.val_; return *this; ,, IN1)
      OPERATOR(-=, val_ -= rhs.val_; return *this; ,, IN1)
      OPERATOR(* , RET(val_ * rhs), const, FIN)
      OPERATOR(/, RET(val_ / rhs), const, FIN)
      OPERATOR(+, RET(val_ + rhs.val_), const, IN1)
      OPERATOR(-, RET(val_ - rhs.val_), const, IN1)
      OPERATOR(==, B_RET(fabsl(val_ - rhs.val_) < 1E-15), const, IN1)
      OPERATOR(!=, B_RET(fabsl(val_ - rhs.val_) >= 1E-15), const, IN1)
      OPERATOR(<, B_RET(val_ < rhs.val_), const, IN1)
      OPERATOR(>, B_RET(val_ > rhs.val_), const, IN1)
      OPERATOR(<=, B_RET(val_ <= rhs.val_), const, IN1)
      OPERATOR(>=, B_RET(val_ >= rhs.val_), const, IN1)
   };
   TMPLT OPERATOR(*, RET(rhs.getSI() * lhs),, const Value lhs, IN1)
   TMPLT OPERATOR(/, INV_RET(lhs / rhs.getSI()),, const Value lhs, IN1)

   template <bool b>
   struct PowCheck
   {
      template <typename T>
      static auto get(const T& ret) { return ret; }
   };

   //Convert to dimensionless numbers if power is not odd
   template <>
   struct PowCheck<false>
   {
      template <typename T>
      static auto get(const T& ret) { return ret.getSI(); }
   };

   TMPLT inline auto sqrt(IN1)
   {
      return PowCheck<!(l % 2 || m % 2 || t % 2 || i % 2
         || k % 2 || n % 2 || j % 2 || a % 2)>::get(SQRT_RET(sqrtl(rhs.getSI())));
   }

   TMPLT inline std::string to_string(IN1)
   {
      return std::string(rhs.getName()) + 
         std::string(": ") +
         std::to_string(rhs.getSI()) + std::string(" ") +
         std::string(rhs.getUnit());
   }

   //Boltzmann Constant
   static const BoltsmanConstant          kb(1.38064852e-23L);
   // Coulomb's constant
   static const ElConst                   e0(8.854187817e-12L);
   // Deal gas constant
   static const GasConstant               R(8.3144598);
   // The gravity of Earth
   static const Acceleration              g(9.80655L);
   // Gravity constant
   static const GrConstant                G(6.674e-11L);
   // Planck constant 
   static const PlanckConstant            h(6.626070040e-34L);
   // Speed of light
   static const Velocity                  c(299792458L);
   // Stefan–Boltzmann constant
   static const StefansConst              s(5.67032e+8L);
   //The Avogadro constant
   static const AvogadroConst             Na(6.02214076e+23);
   // The Faraday electolis constant
   static const FaradayConstant           Fc(96485.33289L);
   // The mass of the electron
   static const Mass                      eMass(9.10938356e-31L);
   // The mass of the neitron
   static const Mass                      nMass(1.674927471e-27L);
   // The charge of the electron
   static const Charge                    eCharge(-1.6021766208e-19L);
   // Rydberg constant
   static const WaveNumber                Rd(10973731.568508L);
   // Vacuum permeability
   static const MgConst                   mu0(1.25663706e-6);

   //SI literals
   LITERALS(Length, m, 1)                 // Metre
   LITERALS(Time, s, 1)                   // Second
   LITERAL(Mass, g, 1e-3)                 // Kilogram
   LITERAL(Mass, kg, 1)                   // Kilogram
   LITERALS(ElCurrent, A, 1)              // Ampere  
   LITERALS(Temperature, K, 1)            // Kelvin   
   LITERALS(LumIntensity, cd, 1)          // candela
   LITERALS(Angle, rad, 1)                // Radian
   LITERALS(Angle, sr, 1)                 // Steradian
   LITERALS(Velocity, m_p_s, 1)           // Metre per second
   LITERALS(Force, N, 1)                  // Newton
   LITERALS(Acceleration, m_p_sqs, 1)     // Metre per square second
   LITERALS(Energy, J, 1)                 // Joule
   LITERALS(Pressure, Pa, 1)              // Pascal
   LITERALS(Power, W, 1)                  // Watt
   LITERALS(ElPotential, V, 1)            // Volt
   LITERALS(ElCapacity, F, 1)             // Farad
   LITERALS(ElInductance, H, 1)           // Henry
   LITERALS(Frequency, Hz, 1)             // Hertz
   LITERALS(Charge, C, 1)                 // Coulomb
   LITERALS(ElResistance, Ohm, 1)         // Ohm
   LITERALS(ElConductance, S, 1)          // Siemens
   LITERALS(Luminance, lx, 1)             // Lux
   LITERALS(DynViscosity, P, 0.1)         // Poise 
   LITERALS(DynViscosity, Pl, 1)          // Poiseuille 
   LITERALS(Energy, eV, 1.6021766208e-19) // Electronvolt
   LITERALS(Energy, Erg, 1e-7)            // Erg
   LITERALS(Energy, cal, 4.1868)          // Calorie
   LITERAL(Time, min, 60)                 // minute
   LITERAL(Length, cm, 1e-2)              // Centimetre
   LITERALS(Pressure, Torr, 133.3)        // Torr
   LITERALS(Pressure, at, 9.80665e+4)     // at
   LITERALS(Pressure, atm, 1.01325e+5)    // atm
   LITERAL(Density, inv_cub_m, 1)         // Inverse cubic metre
   LITERAL(Density, inv_cub_cm, 1e+6)     // Inverse cubic centimetre
   LITERAL(Density, inv_cub_mm, 1e+9)     // Inverse cubic millimetre
   LITERAL(Area, sq_m, 1)                 // Square meter
   LITERAL(Volume, cub_m, 1)              // Cubic meters
   LITERAL(Volume, cub_cm, 1e-6)          // Cubic centimetre
   LITERAL(Volume, cub_mm, 1e-9)          // Cubuc millimetre
   LITERALS(Volume, L, 0.001)             // Litre
   LITERAL(MgStrength, T, 1)              // Tesla
   LITERALS(Activity, Bq, 1)              // Becquerel
   LITERALS(Activity, Ci, 3.7e+10)        // Curie
   LITERALS(Activity, Rd, 1e+6)           // Rutherford
   LITERALS(ChargeMassRatio, R, 2.58e-4)  // Roentgen 
   LITERALS(RadAbsDose, Rad, 0.01)        // Rad
   LITERALS(RadAbsDose, Gy, 1)            // Gray
   LITERALS(RadAbsDose, rem, 0.01)        // Roentgen equivalent man
   LITERALS(RadAbsDose, Sv, 1)            // Sievert
}

#endif //UNITS_H_
