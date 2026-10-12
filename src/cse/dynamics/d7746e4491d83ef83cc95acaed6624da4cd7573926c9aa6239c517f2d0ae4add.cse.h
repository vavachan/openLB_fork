/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-24 Adrian Kummerlaender, Shota Ito
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x21 = cell[15] + cell[17];
auto x22 = cell[12] + x21;
auto x23 = cell[11] + cell[18];
auto x24 = cell[10] + cell[14] + cell[16];
auto x25 = cell[2] + cell[8] + cell[9];
auto x26 = cell[13] + cell[3];
auto x27 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + cell[7] + x22 + x23 + x24 + x25 + x26 + V{1};
auto x28 = V{1} / (x27);
auto x29 = V{0.333333333333333}*cell[13];
auto x30 = V{0.333333333333333}*cell[14];
auto x31 = V{0.333333333333333}*cell[4];
auto x32 = V{0.333333333333333}*cell[5];
auto x33 = V{1}*x28;
auto x34 = -cell[18];
auto x35 = -cell[3];
auto x36 = -cell[8];
auto x37 = cell[9] + x36;
auto x38 = x34 + x35 + x37;
auto x39 = -cell[6];
auto x40 = cell[7] + x39;
auto x41 = -cell[16] + x40;
auto x42 = x22 + x38 + x41;
auto x43 = ((x42)*(x42));
auto x44 = V{0.333333333333333}*cell[0];
auto x45 = V{0.333333333333333}*cell[10];
auto x46 = V{0.333333333333333}*cell[1];
auto x47 = -V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x44 + x45 + x46;
auto x48 = V{0.333333333333333}*cell[11];
auto x49 = V{0.333333333333333}*cell[2];
auto x50 = -V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7] + x48 + x49;
auto x51 = V{0.333333333333333}*cell[15];
auto x52 = V{0.333333333333333}*cell[16];
auto x53 = V{0.333333333333333}*cell[6];
auto x54 = V{0.333333333333333}*cell[7];
auto x55 = cell[13] + cell[17];
auto x56 = -cell[2];
auto x57 = -cell[9];
auto x58 = x23 + x36 + x56 + x57;
auto x59 = -cell[4];
auto x60 = cell[5] + x59;
auto x61 = -cell[14] + x60;
auto x62 = x55 + x58 + x61;
auto x63 = ((x62)*(x62));
auto x64 = V{0.333333333333333}*cell[12];
auto x65 = V{0.333333333333333}*cell[3];
auto x66 = -V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5] + x64 + x65;
auto x67 = V{0.333333333333333}*cell[17];
auto x68 = V{0.333333333333333}*cell[18];
auto x69 = V{0.333333333333333}*cell[8];
auto x70 = V{0.333333333333333}*cell[9];
auto x71 = cell[13] + cell[15];
auto x72 = -cell[1];
auto x73 = -cell[7];
auto x74 = x39 + x72 + x73;
auto x75 = -cell[5] + x59;
auto x76 = x24 + x71 + x74 + x75;
auto x77 = ((x76)*(x76));
auto x78 = x28*x76;
auto x79 = -cell[15] + cell[16];
auto x80 = -cell[17];
auto x81 = cell[18] + x80;
auto x82 = V{1} / (V{3.00000046417339}*util::sqrt(x28*((x20)*(x20))*util::sqrt(((x40 + x42*x78 + x79)*(x40 + x42*x78 + x79)) + ((x28*x42*x62 + x37 + x81)*(x28*x42*x62 + x37 + x81)) + ((-cell[13] + cell[14] + x60 + x62*x78)*(-cell[13] + cell[14] + x60 + x62*x78)) + V{0.5}*((-V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x33*x63 + x47 + x51 + x52 + x53 + x54 + x66)*(-V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x33*x63 + x47 + x51 + x52 + x53 + x54 + x66)) + V{0.5}*((-V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x29 + x30 + x31 + x32 + x33*x43 + x47 + x50)*(-V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x29 + x30 + x31 + x32 + x33*x43 + x47 + x50)) + V{0.5}*((-V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x33*x77 + x44 + x50 + x66 + x67 + x68 + x69 + x70)*(-V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x33*x77 + x44 + x50 + x66 + x67 + x68 + x69 + x70))) + V{0.0277777691819762}/((x19)*(x19))) + V{0.5}/x19);
auto x83 = V{1} / ((x27)*(x27));
auto x84 = V{1.5}*x83;
auto x85 = x77*x84;
auto x86 = x63*x84;
auto x87 = x43*x84;
auto x88 = x86 + x87 + V{-1};
auto x89 = x85 + x88;
auto x90 = V{1} - x82;
auto x91 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x92 = V{3}*cell[14];
auto x93 = V{3}*cell[16];
auto x94 = V{3}*cell[5];
auto x95 = V{3}*cell[7];
auto x96 = V{3}*cell[13] - V{3}*cell[4];
auto x97 = V{3}*cell[15] - V{3}*cell[6];
auto x98 = x28*(V{3}*cell[10] - V{3}*cell[1] + x92 + x93 - x94 - x95 + x96 + x97);
auto x99 = V{3}*x83;
auto x100 = x77*x99;
auto x101 = V{3}*cell[18];
auto x102 = V{3}*cell[9];
auto x103 = V{3}*cell[17] - V{3}*cell[8];
auto x104 = x28*(V{3}*cell[11] - V{3}*cell[2] + x101 - x102 + x103 - x92 + x94 + x96);
auto x105 = x63*x99;
auto x106 = x85 + V{-1};
auto x107 = x28*(V{3}*cell[12] - V{3}*cell[3] - x101 + x102 + x103 - x93 + x95 + x97);
auto x108 = x43*x99;
auto x109 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x110 = V{4.5}*x83;
auto x111 = cell[10] + cell[16] + x74;
auto x112 = x110*((2*cell[13] - 2*cell[4] + x111 + x21 + x58)*(2*cell[13] - 2*cell[4] + x111 + x21 + x58));
auto x113 = x89 + x98;
auto x114 = -x104;
auto x115 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x111 + x25 + x34 + x80;
auto x116 = cell[10] + cell[14] + x72 + x75;
auto x117 = x110*((cell[12] + 2*cell[15] - 2*cell[6] + x116 + x38 + x55)*(cell[12] + 2*cell[15] - 2*cell[6] + x116 + x38 + x55));
auto x118 = -x107;
auto x119 = -cell[12] + x26;
auto x120 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x116 + x119 + x57 + x81;
auto x121 = cell[11] + x56 + x61;
auto x122 = x110*((cell[12] + 2*cell[17] - 2*cell[8] + x121 + x35 + x41 + x71)*(cell[12] + 2*cell[17] - 2*cell[8] + x121 + x35 + x41 + x71));
auto x123 = x104 + x89;
auto x124 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x119 + x121 + x73 + x79;
auto x125 = -x86;
auto x126 = V{1} - x87;
auto x127 = x125 + x126;
auto x128 = x127 + x98;
auto x129 = -x85;
auto x130 = x104 + x129;
auto x131 = x107 + x129;
auto x132 = -x98;
auto x133 = x107 + x89;
auto x0 = V{1}*cell[0]*x90 - x82*(x89*(x29 + x30 + x31 + x32 + x44 + x45 + x46 + x48 + x49 + x51 + x52 + x53 + x54 + x64 + x65 + x67 + x68 + x69 + x70 + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = V{1}*cell[1]*x90 - x82*(x91*(-x100 + x88 + x98) + V{0.0555555555555556});
auto x2 = V{1}*cell[2]*x90 - x82*(x91*(x104 - x105 + x106 + x87) + V{0.0555555555555556});
auto x3 = V{1}*cell[3]*x90 - x82*(x91*(x106 + x107 - x108 + x86) + V{0.0555555555555556});
auto x4 = V{1}*cell[4]*x90 - x82*(x109*(x104 - x112 + x113) + V{0.0277777777777778});
auto x5 = V{1}*cell[5]*x90 - x82*(x109*(-x110*((x115)*(x115)) + x113 + x114) + V{0.0277777777777778});
auto x6 = V{1}*cell[6]*x90 - x82*(x109*(x107 + x113 - x117) + V{0.0277777777777778});
auto x7 = V{1}*cell[7]*x90 - x82*(x109*(-x110*((x120)*(x120)) + x113 + x118) + V{0.0277777777777778});
auto x8 = V{1}*cell[8]*x90 - x82*(x109*(x107 - x122 + x123) + V{0.0277777777777778});
auto x9 = V{1}*cell[9]*x90 - x82*(x109*(-x110*((x124)*(x124)) + x118 + x123) + V{0.0277777777777778});
auto x10 = V{1}*cell[10]*x90 + x82*(x91*(x100 + x128) + V{-0.0555555555555556});
auto x11 = V{1}*cell[11]*x90 + x82*(x91*(x105 + x126 + x130) + V{-0.0555555555555556});
auto x12 = V{1}*cell[12]*x90 + x82*(x91*(x108 + x125 + x131 + V{1}) + V{-0.0555555555555556});
auto x13 = V{1}*cell[13]*x90 + x82*(x109*(x112 + x128 + x130) + V{-0.0277777777777778});
auto x14 = V{1}*cell[14]*x90 - x82*(x109*(-x110*((x115)*(x115)) + x123 + x132) + V{0.0277777777777778});
auto x15 = V{1}*cell[15]*x90 + x82*(x109*(x117 + x128 + x131) + V{-0.0277777777777778});
auto x16 = V{1}*cell[16]*x90 - x82*(x109*(-x110*((x120)*(x120)) + x132 + x133) + V{0.0277777777777778});
auto x17 = V{1}*cell[17]*x90 + x82*(x109*(x107 + x122 + x127 + x130) + V{-0.0277777777777778});
auto x18 = V{1}*cell[18]*x90 - x82*(x109*(-x110*((x124)*(x124)) + x114 + x133) + V{0.0277777777777778});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x27, V{1}*x83*(x43 + x63 + x77) };
}
};

}

}
