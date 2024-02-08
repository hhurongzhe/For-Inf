#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

// This is a code to evaluate interaction matrix elements under plane-wave basis using aPWP method,
// where "potential_chiral_nn" is the primary function.
// In main(), I provide an example to show how to use it for N2LO-opt interaction.

// Constants.
constexpr double gA = 1.29;
constexpr double fpi = 92.4;
constexpr double Mp = 938.2720;
constexpr double Mn = 939.5654;
constexpr double Mnuc = 938.9183;
constexpr double Mpi = 138.0390;
constexpr double Mpi_charged = 139.5702;
constexpr double Mpi_neutral = 134.9766;
constexpr double hbarc = 197.32698;
constexpr double PI = 3.141592653589793;
constexpr double twopicubic = 248.0502134423986;

// LECs structure in operator basis(up to N2LO).
struct LECs_op
{
    double Lambda;
    double Lambda_tilde;
    double C_S_nn;
    double C_T_nn;
    double C_1;
    double C_2;
    double C_3;
    double C_4;
    double C_5;
    double C_6;
    double C_7;
    double D_1;
    double D_2;
    double D_3;
    double D_4;
    double D_5;
    double D_6;
    double D_7;
    double D_8;
    double D_9;
    double D_10;
    double D_11;
    double D_12;
    double D_13;
    double D_14;
    double D_15;
    double c_1;
    double c_3;
    double c_4;
};

// LECs structure in partial-wave basis(up to N2LO).
struct LECs_pw
{
    double Lambda;
    double Lambda_tilde;
    double CC_1s0_nn;
    double CC_3s1_nn;
    double C_1s0;
    double C_3p0;
    double C_1p1;
    double C_3p1;
    double C_3s1;
    double C_3sd1;
    double C_3p2;
    double DD_1s0;
    double D_1s0;
    double D_3p0;
    double D_1p1;
    double D_3p1;
    double DD_3s1;
    double D_3s1;
    double D_3d1;
    double DD_3sd1;
    double D_3sd1;
    double D_1d2;
    double D_3d2;
    double D_3p2;
    double D_3pf2;
    double D_3d3;
    double c_1;
    double c_3;
    double c_4;
};

using lecs_tuple = std::tuple<double, double, double, double, double, double, double, double, double>;

// transformation of LECs between different basis.
lecs_tuple lecs_from_op_to_pw(const lecs_tuple lecs_op)
{
    auto [CS, CT, C1, C2, C3, C4, C5, C6, C7] = lecs_op;
    double fourpi = 4.0 * PI;
    double CC1s0 = fourpi * (CS - 3.0 * CT);
    double CC3s1 = fourpi * (CS + CT);
    double C1s0 = fourpi * (C1 + C2 / 4.0 - 3.0 * C3 - 3.0 / 4.0 * C4 - C6 - C7 / 4.0);
    double C3p0 = fourpi * (-2.0 / 3.0 * C1 + 1.0 / 6.0 * C2 - 2.0 / 3.0 * C3 + C4 / 6.0 - 2.0 / 3.0 * C5 + 2.0 * C6 - C7 / 2.0);
    double C1p1 = fourpi * (-2.0 / 3.0 * C1 + 1.0 / 6.0 * C2 + 2.0 * C3 - C4 / 2.0 + 2.0 / 3.0 * C6 - C7 / 6.0);
    double C3p1 = fourpi * (-2.0 / 3.0 * C1 + 1.0 / 6.0 * C2 - 2.0 / 3.0 * C3 + C4 / 6.0 - 1.0 / 3.0 * C5 - 4.0 / 3.0 * C6 + C7 / 3.0);
    double C3s1 = fourpi * (C1 + C2 / 4.0 + C3 + 1.0 / 4.0 * C4 + C6 / 3.0 + C7 / 12.0);
    double C3sd1 = fourpi * (-2.0 * sqrt(2.0) / 3.0 * C6 - sqrt(2.0) / 6.0 * C7);
    double C3p2 = fourpi * (-2.0 / 3.0 * C1 + C2 / 6.0 - 2.0 / 3.0 * C3 + C4 / 6.0 + C5 / 3.0);
    auto lecs_pw = std::make_tuple(CC1s0, CC3s1, C1s0, C3p0, C1p1, C3p1, C3s1, C3sd1, C3p2);
    return lecs_pw;
}

lecs_tuple lecs_from_pw_to_op(const lecs_tuple lecs_pw)
{
    auto [CC1s0, CC3s1, C1s0, C3p0, C1p1, C3p1, C3s1, C3sd1, C3p2] = lecs_pw;
    double fourpi = 4.0 * PI;
    double CS = (CC1s0 / 4. + (3 * CC3s1) / 4.) / fourpi;
    double CT = (-0.25 * CC1s0 + CC3s1 / 4.) / fourpi;
    double C1 = ((-3 * C1p1) / 16. + C1s0 / 8. - C3p0 / 16. - (3 * C3p1) / 16. -
                 (5 * C3p2) / 16. + (3 * C3s1) / 8.) /
                fourpi;
    double C2 = ((3 * C1p1) / 4. + C1s0 / 2. + C3p0 / 4. + (3 * C3p1) / 4. +
                 (5 * C3p2) / 4. + (3 * C3s1) / 2.) /
                fourpi;
    double C3 = ((3 * C1p1) / 16. - C1s0 / 8. - C3p0 / 16. - C3p2 / 8. + C3s1 / 8. +
                 C3sd1 / (4. * sqrt(2))) /
                fourpi;
    double C4 = ((-3 * C1p1) / 4. - C1s0 / 2. + C3p0 / 4. + C3p2 / 2. + C3s1 / 2. +
                 C3sd1 / sqrt(2)) /
                fourpi;
    double C5 = (-0.5 * C3p0 - (3 * C3p1) / 4. + (5 * C3p2) / 4.) / fourpi;
    double C6 = (C3p0 / 8. - (3 * C3p1) / 16. + C3p2 / 16. -
                 (3 * C3sd1) / (4. * sqrt(2))) /
                fourpi;
    double C7 = (-0.5 * C3p0 + (3 * C3p1) / 4. - C3p2 / 4. - (3 * C3sd1) / sqrt(2)) /
                fourpi;
    auto lecs_op = std::make_tuple(CS, CT, C1, C2, C3, C4, C5, C6, C7);
    return lecs_op;
}

// automated plane-wave projection.
double potential_auto(double f1, double f2, double f3, double f4, double f5, double f6,
                      double ppx, double ppy, double ppz, double px, double py, double pz,
                      int s1_out, int s2_out, int s1_in, int s2_in)
{
    // params f_j : coefficients of operators w_j
    // params ppx, ppy, ppz : final relative momentum p', in MeV
    // params px, py, pz : initial relative momentum p, in MeV
    // params s1_out, s2_out : final spin of 2 nucleons(doubled)
    // params s1_in, s2_in : initial spin of 2 nucleons(doubled)
    // return : interaction matrix element, in MeV^(-2)

    std::tuple<int, int, int, int> spin_tuple(s1_out, s2_out, s1_in, s2_in);

    if (spin_tuple == std::make_tuple(1, 1, 1, 1))
    {
        return f1 + f2 + f5 * pow(ppz, 2) + f6 * pow(ppz, 2) + f4 * pow(ppy, 2) * pow(px, 2) -
               2 * f4 * ppx * ppy * px * py + f4 * pow(ppx, 2) * pow(py, 2) + 2 * f5 * ppz * pz - 2 * f6 * ppz * pz +
               f5 * pow(pz, 2) + f6 * pow(pz, 2);
    }
    else if (spin_tuple == std::make_tuple(1, 1, 1, -1))
    {
        return -(f3 * ppz * px) + f4 * ppy * ppz * px * py - f4 * ppx * ppz * pow(py, 2) + f6 * (ppx - px) * (ppz - pz) +
               f3 * ppx * pz - f4 * pow(ppy, 2) * px * pz + f4 * ppx * ppy * py * pz + f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(1, 1, -1, 1))
    {
        return -(f3 * ppz * px) + f4 * ppy * ppz * px * py - f4 * ppx * ppz * pow(py, 2) + f6 * (ppx - px) * (ppz - pz) +
               f3 * ppx * pz - f4 * pow(ppy, 2) * px * pz + f4 * ppx * ppy * py * pz + f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(1, 1, -1, -1))
    {
        return f5 * (pow(ppx, 2) - pow(ppy, 2) + 2 * ppx * px + pow(px, 2) - 2 * ppy * py - pow(py, 2)) +
               f6 * (pow(ppx, 2) - pow(ppy, 2) - 2 * ppx * px + pow(px, 2) + 2 * ppy * py -
                     pow(py, 2)) +
               f4 * (pow(ppz, 2) * (-pow(px, 2) + pow(py, 2)) +
                     2 * ppz * (ppx * px - ppy * py) * pz + (-pow(ppx, 2) + pow(ppy, 2)) * pow(pz, 2));
    }
    else if (spin_tuple == std::make_tuple(1, -1, 1, 1))
    {
        return f3 * ppz * px + f4 * ppy * ppz * px * py - f4 * ppx * ppz * pow(py, 2) + f6 * (ppx - px) * (ppz - pz) -
               f3 * ppx * pz - f4 * pow(ppy, 2) * px * pz + f4 * ppx * ppy * py * pz + f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(1, -1, 1, -1))
    {
        return f1 - f2 - f5 * pow(ppz, 2) - f6 * pow(ppz, 2) - f4 * pow(ppy, 2) * pow(px, 2) +
               2 * f4 * ppx * ppy * px * py - f4 * pow(ppx, 2) * pow(py, 2) - 2 * f5 * ppz * pz + 2 * f6 * ppz * pz -
               f5 * pow(pz, 2) - f6 * pow(pz, 2);
    }
    else if (spin_tuple == std::make_tuple(1, -1, -1, 1))
    {
        return 2 * f2 + f6 * pow(ppx, 2) + f6 * pow(ppy, 2) - 2 * f6 * ppx * px + f6 * pow(px, 2) +
               f4 * pow(ppz, 2) * pow(px, 2) - 2 * f6 * ppy * py + f6 * pow(py, 2) +
               f4 * pow(ppz, 2) * pow(py, 2) + f5 * (pow(ppx, 2) + pow(ppy, 2) + 2 * ppx * px + pow(px, 2) + 2 * ppy * py + pow(py, 2)) -
               2 * f4 * ppx * ppz * px * pz - 2 * f4 * ppy * ppz * py * pz + f4 * pow(ppx, 2) * pow(pz, 2) +
               f4 * pow(ppy, 2) * pow(pz, 2);
    }
    else if (spin_tuple == std::make_tuple(1, -1, -1, -1))
    {
        return -(f3 * ppz * px) - f4 * ppy * ppz * px * py + f4 * ppx * ppz * pow(py, 2) - f6 * (ppx - px) * (ppz - pz) +
               f3 * ppx * pz + f4 * pow(ppy, 2) * px * pz - f4 * ppx * ppy * py * pz - f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(-1, 1, 1, 1))
    {
        return f3 * ppz * px + f4 * ppy * ppz * px * py - f4 * ppx * ppz * pow(py, 2) + f6 * (ppx - px) * (ppz - pz) -
               f3 * ppx * pz - f4 * pow(ppy, 2) * px * pz + f4 * ppx * ppy * py * pz + f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(-1, 1, 1, -1))
    {
        return 2 * f2 + f6 * pow(ppx, 2) + f6 * pow(ppy, 2) - 2 * f6 * ppx * px + f6 * pow(px, 2) +
               f4 * pow(ppz, 2) * pow(px, 2) - 2 * f6 * ppy * py + f6 * pow(py, 2) +
               f4 * pow(ppz, 2) * pow(py, 2) + f5 * (pow(ppx, 2) + pow(ppy, 2) + 2 * ppx * px + pow(px, 2) + 2 * ppy * py + pow(py, 2)) -
               2 * f4 * ppx * ppz * px * pz - 2 * f4 * ppy * ppz * py * pz + f4 * pow(ppx, 2) * pow(pz, 2) +
               f4 * pow(ppy, 2) * pow(pz, 2);
    }
    else if (spin_tuple == std::make_tuple(-1, 1, -1, 1))
    {
        return f1 - f2 - f5 * pow(ppz, 2) - f6 * pow(ppz, 2) - f4 * pow(ppy, 2) * pow(px, 2) +
               2 * f4 * ppx * ppy * px * py - f4 * pow(ppx, 2) * pow(py, 2) - 2 * f5 * ppz * pz + 2 * f6 * ppz * pz -
               f5 * pow(pz, 2) - f6 * pow(pz, 2);
    }
    else if (spin_tuple == std::make_tuple(-1, 1, -1, -1))
    {
        return -(f3 * ppz * px) - f4 * ppy * ppz * px * py + f4 * ppx * ppz * pow(py, 2) - f6 * (ppx - px) * (ppz - pz) +
               f3 * ppx * pz + f4 * pow(ppy, 2) * px * pz - f4 * ppx * ppy * py * pz - f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(-1, -1, 1, 1))
    {
        return f5 * (pow(ppx, 2) - pow(ppy, 2) + 2 * ppx * px + pow(px, 2) - 2 * ppy * py - pow(py, 2)) +
               f6 * (pow(ppx, 2) - pow(ppy, 2) - 2 * ppx * px + pow(px, 2) + 2 * ppy * py -
                     pow(py, 2)) +
               f4 * (pow(ppz, 2) * (-pow(px, 2) + pow(py, 2)) +
                     2 * ppz * (ppx * px - ppy * py) * pz + (-pow(ppx, 2) + pow(ppy, 2)) * pow(pz, 2));
    }
    else if (spin_tuple == std::make_tuple(-1, -1, 1, -1))
    {
        return f3 * ppz * px - f4 * ppy * ppz * px * py + f4 * ppx * ppz * pow(py, 2) - f6 * (ppx - px) * (ppz - pz) -
               f3 * ppx * pz + f4 * pow(ppy, 2) * px * pz - f4 * ppx * ppy * py * pz - f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(-1, -1, -1, 1))
    {
        return f3 * ppz * px - f4 * ppy * ppz * px * py + f4 * ppx * ppz * pow(py, 2) - f6 * (ppx - px) * (ppz - pz) -
               f3 * ppx * pz + f4 * pow(ppy, 2) * px * pz - f4 * ppx * ppy * py * pz - f5 * (ppx + px) * (ppz + pz);
    }
    else if (spin_tuple == std::make_tuple(-1, -1, -1, -1))
    {
        return f1 + f2 + f5 * pow(ppz, 2) + f6 * pow(ppz, 2) + f4 * pow(ppy, 2) * pow(px, 2) -
               2 * f4 * ppx * ppy * px * py + f4 * pow(ppx, 2) * pow(py, 2) + 2 * f5 * ppz * pz - 2 * f6 * ppz * pz +
               f5 * pow(pz, 2) + f6 * pow(pz, 2);
    }
    else
    {
        throw std::invalid_argument("Invalid spin configuration in function potential_auto!");
    }
}

// LO contact terms.
std::vector<double> potential_lo_contact(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);

    f[0] = f[0] + lecs.C_S_nn;
    f[1] = f[1] + lecs.C_T_nn;

    double Lambda = lecs.Lambda;
    int n = 3;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// NLO contact terms.
std::vector<double> potential_nlo_contact(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double kx = (ppx + px) / 2;
    double ky = (ppy + py) / 2;
    double kz = (ppz + pz) / 2;
    double q2 = qx * qx + qy * qy + qz * qz;
    double k2 = kx * kx + ky * ky + kz * kz;

    f[0] = f[0] + lecs.C_1 * q2 + lecs.C_2 * k2;
    f[1] = f[1] + lecs.C_3 * q2 + lecs.C_4 * k2;
    f[2] = f[2] + lecs.C_5 / 2.0;
    f[4] = f[4] + lecs.C_7 / 4.0;
    f[5] = f[5] + lecs.C_6;

    double Lambda = lecs.Lambda;
    int n = 2;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// N3LO contact terms, not used by now.
std::vector<double> potential_n3lo_contact(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double kx = (ppx + px) / 2.0;
    double ky = (ppy + py) / 2.0;
    double kz = (ppz + pz) / 2.0;
    double q2 = qx * qx + qy * qy + qz * qz;
    double k2 = kx * kx + ky * ky + kz * kz;
    double q4 = q2 * q2;
    double k4 = k2 * k2;
    double qk2 = pow(ky * qx - kx * qy, 2) + pow(-(kz * qx) + kx * qz, 2) +
                 pow(kz * qy - ky * qz, 2);
    f[0] = f[0] + lecs.D_1 * q4 + lecs.D_2 * k4 + lecs.D_3 * q2 * k2 + lecs.D_4 * qk2;
    f[1] = f[1] + lecs.D_5 * q4 + lecs.D_6 * k4 + lecs.D_7 * q2 * k2 + lecs.D_8 * qk2;
    f[2] = f[2] + (lecs.D_9 * q2 + lecs.D_10 * k2) / 2.0;
    f[3] = f[3] + lecs.D_15;
    f[4] = f[4] + (lecs.D_13 * q2 + lecs.D_14 * k2) / 4.0;
    f[5] = f[5] + lecs.D_11 * q2 + lecs.D_12 * k2;

    double Lambda = lecs.Lambda;
    int n = 2;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// one-neutral-pion exchange term.
std::vector<double> potential_one_pion_exchange_neutral(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double q2 = qx * qx + qy * qy + qz * qz;
    double prefactor = -1.0 * gA * gA / (4 * fpi * fpi);
    f[5] = f[5] + prefactor / (q2 + Mpi_neutral * Mpi_neutral);

    double Lambda = lecs.Lambda;
    int n = 4;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// one-charged-pion exchange term.
std::vector<double> potential_one_pion_exchange_charged(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double q2 = qx * qx + qy * qy + qz * qz;
    double prefactor = -1.0 * gA * gA / (4 * fpi * fpi);
    f[5] = f[5] + prefactor / (q2 + Mpi_charged * Mpi_charged);

    double Lambda = lecs.Lambda;
    int n = 4;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// Loop function L
double loop_function_L(double qmag, const LECs_op &lecs)
{
    const double eps = 1e-2;
    if (qmag < eps)
    {
        qmag = 1e-2;
    }
    double temp;
    double lambda = lecs.Lambda_tilde;
    double w = sqrt(4.0 * Mpi * Mpi + qmag * qmag);
    double num = lambda * lambda * (2.0 * Mpi * Mpi + qmag * qmag) - 2.0 * Mpi * Mpi * qmag * qmag +
                 lambda * sqrt(lambda * lambda - 4.0 * Mpi * Mpi) * qmag * w;
    double den = 2.0 * Mpi * Mpi * (lambda * lambda + qmag * qmag);
    double fac = w / (2.0 * qmag);
    temp = fac * log(num / den);
    return temp;
}

// Loop function A
double loop_function_A(double qmag, const LECs_op &lecs)
{
    const double eps = 1e-2;
    if (qmag < eps)
    {
        qmag = 1e-2;
    }
    double temp;
    double lambda = lecs.Lambda_tilde;
    double num = qmag * (lambda - 2.0 * Mpi);
    double den = qmag * qmag + 2.0 * lambda * Mpi;
    double fac = 1.0 / (2.0 * qmag);
    temp = fac * atan(num / den);
    return temp;
}

// NLO two-pion exchange term
std::vector<double> potential_nlo_two_pion_exchange(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double q2 = qx * qx + qy * qy + qz * qz;
    double qmag = sqrt(q2);

    double isospin_factor = 1.0; // todo: this factor accounts for tau_1.tau_2 term, which is dependent with total isospin I.
    double f1_fac = loop_function_L(qmag, lecs) / (384.0 * PI * PI * pow(fpi, 4));
    double f1_part1 = 4.0 * pow(Mpi, 2) * (1.0 + 4.0 * pow(gA, 2) - 5.0 * pow(gA, 4));
    double f1_part2 = qmag * qmag * (1.0 + 10.0 * pow(gA, 2) - 23.0 * pow(gA, 4));
    double f1_part3 = -48.0 * pow(gA, 4) * pow(Mpi, 4) / (4.0 * pow(Mpi, 2) + qmag * qmag);
    double f1 = isospin_factor * f1_fac * (f1_part1 + f1_part2 + f1_part3);
    double f6 = -3.0 * pow(gA, 4) / (64.0 * PI * PI * pow(fpi, 4)) * loop_function_L(qmag, lecs);
    double f2 = -qmag * qmag * f6;
    f[0] = f[0] + f1;
    f[1] = f[1] + f2;
    f[5] = f[5] + f6;

    double Lambda = lecs.Lambda;
    int n = 2;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// N2LO two-pion exchange term
std::vector<double> potential_n2lo_two_pion_exchange(double ppx, double ppy, double ppz, double px, double py, double pz, const LECs_op &lecs)
{
    std::vector<double> f(6, 0.0);
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);
    double qx = ppx - px;
    double qy = ppy - py;
    double qz = ppz - pz;
    double q2 = qx * qx + qy * qy + qz * qz;
    double qmag = sqrt(q2);

    double isospin_factor = 1.0; // todo: this factor accounts for tau_1.tau_2 term, which is dependent with total isospin I.
    double gaga = pow(gA, 2);
    double ffff = pow(fpi, 4);
    double mm = pow(Mpi, 2);
    double f1_part1 = 3.0 * gaga / (16.0 * PI * ffff);
    double f1_part2 = 2.0 * mm * (lecs.c_3 - 2.0 * lecs.c_1) + lecs.C_3 * qmag * qmag;
    double f1_part3 = 2.0 * mm + qmag * qmag;
    double f1_part4 = loop_function_A(qmag, lecs);
    double f1 = f1_part1 * f1_part2 * f1_part3 * f1_part4;
    double f6 = -isospin_factor * gaga / (32.0 * PI * ffff) * lecs.c_4 * (4.0 * mm + qmag * qmag) *
                loop_function_A(qmag, lecs);
    double f2 = -qmag * qmag * f6;
    f[0] = f[0] + f1;
    f[1] = f[1] + f2;
    f[5] = f[5] + f6;

    double Lambda = lecs.Lambda;
    int n = 2;
    double regulator = std::exp(-(std::pow(pmag, 2 * n) + std::pow(ppmag, 2 * n)) / std::pow(Lambda, 2 * n));

    for (auto &component : f)
    {
        component *= regulator;
    }

    return f;
}

// total chiral nn potential.
double potential_chiral_nn(const std::vector<double> &p_out, int s1_out, int s2_out,
                           const std::vector<double> &p_in, int s1_in, int s2_in, const LECs_op &lecs)
{
    // params p_out : final relative momentum p' 3-dim vector, in MeV
    // params p_in : initial relative momentum p 3-dim vector, in MeV
    // params s1_out, s2_out : final spin of 2 nucleons(doubled)
    // params s1_in, s2_in : initial spin of 2 nucleons(doubled)
    // params lecs: LECs structure in operator basis
    // return : interaction matrix elements in MeV^(-2)
    double ppx = p_out[0];
    double ppy = p_out[1];
    double ppz = p_out[2];
    double px = p_in[0];
    double py = p_in[1];
    double pz = p_in[2];
    double ppmag = std::sqrt(ppx * ppx + ppy * ppy + ppz * ppz);
    double pmag = std::sqrt(px * px + py * py + pz * pz);

    std::vector<double> component_vec(6, 0.0);

    // lo terms.
    std::vector<double> lo_contact = potential_lo_contact(ppx, ppy, ppz, px, py, pz, lecs);
    std::vector<double> one_pion_exchange = potential_one_pion_exchange_neutral(ppx, ppy, ppz, px, py, pz, lecs);

    // nlo terms.
    std::vector<double> nlo_contact = potential_nlo_contact(ppx, ppy, ppz, px, py, pz, lecs);
    std::vector<double> nlo_two_pion_exchange = potential_nlo_two_pion_exchange(ppx, ppy, ppz, px, py, pz, lecs);

    // n2lo terms.
    std::vector<double> n2lo_two_pion_exchange = potential_n2lo_two_pion_exchange(ppx, ppy, ppz, px, py, pz, lecs);

    for (size_t i = 0; i < component_vec.size(); ++i)
    {
        component_vec[i] += lo_contact[i] + one_pion_exchange[i];
        component_vec[i] += nlo_contact[i] + nlo_two_pion_exchange[i];
        component_vec[i] += n2lo_two_pion_exchange[i];
    }

    // total coefficients.
    double f1 = component_vec[0];
    double f2 = component_vec[1];
    double f3 = component_vec[2];
    double f4 = component_vec[3];
    double f5 = component_vec[4];
    double f6 = component_vec[5];

    // minimal-relativity term.
    double e = std::sqrt(Mnuc * Mnuc + pmag * pmag);
    double ep = std::sqrt(Mnuc * Mnuc + ppmag * ppmag);
    double min_rel = std::sqrt(Mnuc / e) * std::sqrt(Mnuc / ep);

    return potential_auto(f1, f2, f3, f4, f5, f6, ppx, ppy, ppz, px, py, pz, s1_out, s2_out, s1_in, s2_in) * min_rel / twopicubic;
}

int main()
{
    // N2LO-opt LECs.
    constexpr double Lambda = 500;
    constexpr double Lambda_tilde = 700;
    constexpr double c1 = -0.91863953 * 1e-3;
    constexpr double c3 = -3.88868749 * 1e-3;
    constexpr double c4 = 4.31032716 * 1e-3;
    constexpr double CC1s0nn = -0.15176475 * 1e-2;
    constexpr double CC3s1 = -0.15843418 * 1e-2;
    constexpr double C1s0 = 2.40402194 * 1e-8;
    constexpr double C3p0 = 1.26339076 * 1e-8;
    constexpr double C1p1 = 0.41704554 * 1e-8;
    constexpr double C3p1 = -0.78265850 * 1e-8;
    constexpr double C3s1 = 0.92838466 * 1e-8;
    constexpr double C3sd1 = 0.61814142 * 1e-8;
    constexpr double C3p2 = -0.67780851 * 1e-8;

    // transform to operator basis before aPWP.
    auto lecs_tuple_pw = std::make_tuple(CC1s0nn, CC3s1, C1s0, C3p0, C1p1, C3p1, C3s1, C3sd1, C3p2);
    auto lecs_tuple_op = lecs_from_pw_to_op(lecs_tuple_pw);

    // print LECs.
    std::cout << "operator LECs: " << std::get<0>(lecs_tuple_op) << ", "
              << std::get<1>(lecs_tuple_op) << ", "
              << std::get<2>(lecs_tuple_op) << ", "
              << std::get<3>(lecs_tuple_op) << ", "
              << std::get<4>(lecs_tuple_op) << ", "
              << std::get<5>(lecs_tuple_op) << ", "
              << std::get<6>(lecs_tuple_op) << ", "
              << std::get<7>(lecs_tuple_op) << ", "
              << std::get<8>(lecs_tuple_op) << std::endl;

    // store LECs in the struct.
    LECs_op lecs_op;
    lecs_op.Lambda = Lambda;
    lecs_op.Lambda_tilde = Lambda_tilde;
    lecs_op.c_1 = c1;
    lecs_op.c_3 = c3;
    lecs_op.c_4 = c4;
    lecs_op.C_S_nn = std::get<0>(lecs_tuple_op);
    lecs_op.C_T_nn = std::get<1>(lecs_tuple_op);
    lecs_op.C_1 = std::get<2>(lecs_tuple_op);
    lecs_op.C_2 = std::get<3>(lecs_tuple_op);
    lecs_op.C_3 = std::get<4>(lecs_tuple_op);
    lecs_op.C_4 = std::get<5>(lecs_tuple_op);
    lecs_op.C_5 = std::get<6>(lecs_tuple_op);
    lecs_op.C_6 = std::get<7>(lecs_tuple_op);
    lecs_op.C_7 = std::get<8>(lecs_tuple_op);

    // 3-dim momentum vectors p, p'.
    std::vector<double> p_out = {0, 200, 300};
    std::vector<double> p_in = {100, 200, 300};

    for (int s1_out = -1; s1_out <= 1; s1_out = s1_out + 2)
    {
        for (int s2_out = -1; s2_out <= 1; s2_out = s2_out + 2)
        {
            for (int s1_in = -1; s1_in <= 1; s1_in = s1_in + 2)
            {
                for (int s2_in = -1; s2_in <= 1; s2_in = s2_in + 2)
                {
                    double mtx = potential_chiral_nn(p_out, s1_out, s2_out, p_in, s1_in, s2_in, lecs_op);
                    std::cout << s1_out << "," << s2_out << "," << s1_in << "," << s2_in << "," << mtx << std::endl;
                }
            }
        }
    }

    return 0;
}
