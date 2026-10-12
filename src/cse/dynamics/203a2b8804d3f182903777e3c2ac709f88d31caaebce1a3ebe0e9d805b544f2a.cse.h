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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::P1Momentum, momenta::ZeroStress, momenta::DefineToNEq>, equilibria::P1, collision::P1, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x28 = parameters.template get<collision::P1::ABSORPTION>();
auto x27 = parameters.template get<collision::P1::SCATTERING>();
auto x29 = V{0.0952380952380952}*cell[7];
auto x30 = V{0.0476190476190476}*cell[0];
auto x31 = V{0.0476190476190476}*cell[16];
auto x32 = V{0.0476190476190476}*cell[3];
auto x33 = V{0.0952380952380952}*cell[10];
auto x34 = V{0.0952380952380952}*cell[11];
auto x35 = V{0.0952380952380952}*cell[4] + x30 + x31 + x32 + x33 + x34;
auto x36 = V{0.0476190476190476}*cell[15];
auto x37 = V{0.0476190476190476}*cell[2];
auto x38 = V{0.0952380952380952}*cell[12];
auto x39 = V{0.0952380952380952}*cell[6] + x36 + x37 + x38;
auto x40 = V{0.0952380952380952}*cell[13];
auto x41 = V{0.0952380952380952}*cell[5] + x40;
auto x42 = V{0.0476190476190476}*cell[21] + V{0.0476190476190476}*cell[22] + V{0.0476190476190476}*cell[8] + V{0.0476190476190476}*cell[9];
auto x43 = V{1}*x27;
auto x44 = V{1}*x28;
auto x45 = V{0.0952380952380952}*cell[9];
auto x46 = V{0.0952380952380952}*cell[25];
auto x47 = V{0.0952380952380952}*cell[26];
auto x48 = V{0.0952380952380952}*cell[18] + x46 + x47;
auto x49 = V{0.0476190476190476}*cell[14];
auto x50 = V{0.0476190476190476}*cell[1];
auto x51 = V{0.0952380952380952}*cell[8] + x49 + x50;
auto x52 = V{0.0476190476190476}*cell[19] + V{0.0476190476190476}*cell[20] + V{0.0476190476190476}*cell[6] + V{0.0476190476190476}*cell[7];
auto x53 = V{0.0952380952380952}*cell[22];
auto x54 = V{0.0952380952380952}*cell[24];
auto x55 = V{0.0952380952380952}*cell[20] + x30 + x54;
auto x56 = V{0.0476190476190476}*cell[17] + V{0.0476190476190476}*cell[18] + V{0.0476190476190476}*cell[4] + V{0.0476190476190476}*cell[5];
auto x57 = V{0.0380952380952381}*cell[17];
auto x58 = V{0.0380952380952381}*cell[0];
auto x59 = V{0.0380952380952381}*cell[13];
auto x60 = V{0.0380952380952381}*cell[26];
auto x61 = V{0.0761904761904762}*cell[9];
auto x62 = x58 + x59 + x60 + x61;
auto x63 = V{0.0380952380952381}*cell[12];
auto x64 = V{0.0380952380952381}*cell[25];
auto x65 = V{0.0761904761904762}*cell[8];
auto x66 = x63 + x64 + x65;
auto x67 = V{0.0380952380952381}*cell[16];
auto x68 = V{0.0380952380952381}*cell[3];
auto x69 = V{0.0761904761904762}*cell[1];
auto x70 = V{0.0761904761904762}*cell[6];
auto x71 = V{0.0761904761904762}*cell[7];
auto x72 = x67 + x68 + x69 + x70 + x71;
auto x73 = V{0.0380952380952381}*cell[24];
auto x74 = V{0.114285714285714}*cell[11] - x73;
auto x75 = V{0.0761904761904762}*cell[2];
auto x76 = V{0.0380952380952381}*cell[23];
auto x77 = V{0.114285714285714}*cell[10] - x76;
auto x78 = x75 + x77;
auto x79 = V{0.0380952380952381}*cell[18];
auto x80 = V{0.0380952380952381}*cell[5];
auto x81 = x79 + x80;
auto x82 = V{1.414214}*x27;
auto x83 = V{1.414214}*x28;
auto x84 = V{0.0761904761904762}*cell[15];
auto x85 = V{0.0761904761904762}*cell[21];
auto x86 = V{0.0761904761904762}*cell[22];
auto x87 = x84 + x85 + x86;
auto x88 = V{0.114285714285714}*cell[12];
auto x89 = -x64;
auto x90 = V{0.0380952380952381}*cell[11];
auto x91 = x58 + x73 + x90;
auto x92 = x88 + x89 + x91;
auto x93 = V{0.114285714285714}*cell[13];
auto x94 = -x60;
auto x95 = V{0.0380952380952381}*cell[10];
auto x96 = x76 + x95;
auto x97 = x93 + x94 + x96;
auto x98 = V{0.0380952380952381}*cell[4];
auto x99 = x57 + x98;
auto x100 = V{0.0380952380952381}*cell[19];
auto x101 = V{0.0761904761904762}*cell[3];
auto x102 = V{0.0761904761904762}*cell[4];
auto x103 = x101 + x102;
auto x104 = V{0.0380952380952381}*cell[15];
auto x105 = V{0.0380952380952381}*cell[20];
auto x106 = V{0.0380952380952381}*cell[2];
auto x107 = V{0.0380952380952381}*cell[7];
auto x108 = x104 + x105 + x106 + x107;
auto x109 = V{0.0761904761904762}*cell[5];
auto x110 = x59 + x60;
auto x111 = x109 + x110;
auto x112 = V{0.0761904761904762}*cell[16];
auto x113 = x112 + x85;
auto x114 = x63 + x64;
auto x115 = x114 + x58;
auto x116 = x102 + x74;
auto x117 = V{0.0380952380952381}*cell[6];
auto x118 = x100 + x104 + x106 + x117;
auto x119 = V{0.0380952380952381}*cell[21];
auto x120 = V{0.0380952380952381}*cell[14];
auto x121 = V{0.0380952380952381}*cell[1];
auto x122 = V{0.0761904761904762}*cell[18];
auto x123 = x120 + x121 + x122;
auto x124 = V{0.0761904761904762}*cell[20];
auto x125 = V{0.114285714285714}*cell[26] - x59;
auto x126 = x124 + x125 + x91;
auto x127 = V{0.0380952380952381}*cell[22];
auto x128 = V{0.0380952380952381}*cell[9];
auto x129 = x114 + x127 + x128;
auto x130 = V{0.0761904761904762}*cell[19];
auto x131 = x110 + x130;
auto x132 = x112 + x71;
auto x133 = V{0.0380952380952381}*cell[8];
auto x134 = x58 + x96;
auto x135 = x119 + x133 + x134;
auto x136 = V{0.114285714285714}*cell[25] - x63;
auto x137 = x136 + x75;
auto x138 = V{0.0642857142857143}*cell[23];
auto x139 = V{0.0321428571428571}*cell[0];
auto x140 = V{0.0321428571428571}*cell[22];
auto x141 = V{0.0321428571428571}*cell[9];
auto x142 = V{0.0642857142857143}*cell[11];
auto x143 = V{0.0642857142857143}*cell[12];
auto x144 = V{0.0642857142857143}*cell[1];
auto x145 = x139 + x140 + x141 + x142 + x143 + x144;
auto x146 = V{0.0321428571428571}*cell[18];
auto x147 = V{0.0321428571428571}*cell[5];
auto x148 = V{0.0642857142857143}*cell[26];
auto x149 = V{0.0642857142857143}*cell[3];
auto x150 = x146 + x147 + x148 + x149;
auto x151 = V{0.0321428571428571}*cell[20];
auto x152 = V{0.0321428571428571}*cell[7];
auto x153 = V{0.0642857142857143}*cell[2];
auto x154 = x151 + x152 + x153;
auto x155 = V{0.0321428571428571}*cell[17];
auto x156 = V{0.0964285714285714}*cell[4] - x155;
auto x157 = V{0.0321428571428571}*cell[19];
auto x158 = V{0.0964285714285714}*cell[6] - x157;
auto x159 = V{0.0321428571428571}*cell[21];
auto x160 = V{0.0964285714285714}*cell[8] - x159;
auto x161 = V{1.732051}*x27;
auto x162 = V{1.732051}*x28;
auto x163 = V{0.0642857142857143}*cell[24];
auto x164 = V{0.0321428571428571}*cell[6];
auto x165 = V{0.0642857142857143}*cell[16];
auto x166 = V{0.0964285714285714}*cell[7] - x151 + x157 + x164 + x165;
auto x167 = V{0.0642857142857143}*cell[10];
auto x168 = V{0.0642857142857143}*cell[13];
auto x169 = V{0.0321428571428571}*cell[8];
auto x170 = x139 + x159 + x169;
auto x171 = x144 + x167 + x168 + x170;
auto x172 = V{0.0642857142857143}*cell[25];
auto x173 = x146 + x147 + x172;
auto x174 = V{0.0964285714285714}*cell[9] - x140;
auto x175 = V{0.0642857142857143}*cell[15];
auto x176 = V{0.0964285714285714}*cell[22] - x141 + x175;
auto x177 = x151 + x152 + x163;
auto x178 = V{0.0321428571428571}*cell[4];
auto x179 = x155 + x178;
auto x180 = x149 + x179;
auto x181 = V{0.0964285714285714}*cell[5] - x146;
auto x182 = x138 + x179;
auto x183 = V{0.0964285714285714}*cell[21] - x169 + x175;
auto x184 = V{0.0952380952380952}*cell[23];
auto x185 = V{0.0952380952380952}*cell[17] + x184 + x31 + x32;
auto x186 = V{0.0952380952380952}*cell[19] + x36 + x37;
auto x187 = V{0.0952380952380952}*cell[21] + x30 + x49 + x50;
auto x188 = V{0.114285714285714}*cell[23] - x95;
auto x189 = V{0.0761904761904762}*cell[14];
auto x190 = x189 + x67 + x68;
auto x191 = V{0.114285714285714}*cell[24] - x90;
auto x192 = x124 + x191;
auto x193 = V{0.0761904761904762}*cell[17];
auto x194 = x188 + x193;
auto x195 = x122 + x189;
auto x196 = x101 + x193;
auto x197 = x120 + x121 + x84;
auto x198 = V{0.0642857142857143}*cell[14];
auto x199 = V{0.0964285714285714}*cell[17] - x178 + x198;
auto x200 = V{0.0964285714285714}*cell[19] - x164 + x165;
auto x201 = x139 + x140 + x141;
auto x202 = V{0.0964285714285714}*cell[20] - x152 + x157 + x164;
auto x203 = V{0.0964285714285714}*cell[18] - x147 + x198;
auto x0 = cell[0];
auto x1 = cell[1] + x43*(-V{0.904761904761905}*cell[1] + x29 + x35 + x39 + x41 + x42) - x44*(cell[1] + V{0.0476190476190476});
auto x2 = cell[2] + x43*(-V{0.904761904761905}*cell[2] + x35 + x45 + x48 + x51 + x52) - x44*(cell[2] + V{0.0476190476190476});
auto x3 = cell[3] + x43*(-V{0.904761904761905}*cell[3] + x33 + x39 + x47 + x51 + x53 + x55 + x56) - x44*(cell[3] + V{0.0476190476190476});
auto x4 = cell[4] + x82*(-V{0.885714285714286}*cell[4] - x57 + x62 + x66 + x72 + x74 + x78 + x81) - x83*(cell[4] + V{0.0380952380952381});
auto x5 = cell[5] + x82*(-V{0.885714285714286}*cell[5] + x72 - x79 + x87 + x92 + x97 + x99) - x83*(cell[5] + V{0.0380952380952381});
auto x6 = cell[6] + x82*(-V{0.885714285714286}*cell[6] - x100 + x103 + x108 + x111 + x65 + x69 + x77 + x86 + x92) - x83*(cell[6] + V{0.0380952380952381});
auto x7 = cell[7] + x82*(-V{0.885714285714286}*cell[7] - x105 + x109 + x113 + x115 + x116 + x118 + x61 + x69 + x97) - x83*(cell[7] + V{0.0380952380952381});
auto x8 = cell[8] + x82*(-V{0.885714285714286}*cell[8] + x103 - x119 + x123 + x126 + x129 + x70 + x78) - x83*(cell[8] + V{0.0380952380952381});
auto x9 = cell[9] + x82*(-V{0.885714285714286}*cell[9] + x116 + x123 - x127 + x131 + x132 + x135 + x137) - x83*(cell[9] + V{0.0380952380952381});
auto x10 = cell[10] + x161*(-V{0.871428571428571}*cell[10] - x138 + x145 + x150 + x154 + x156 + x158 + x160) - x162*(cell[10] + V{0.0321428571428571});
auto x11 = cell[11] + x161*(-V{0.871428571428571}*cell[11] + x153 + x156 - x163 + x166 + x171 + x173 + x174) - x162*(cell[11] + V{0.0321428571428571});
auto x12 = cell[12] + x161*(-V{0.871428571428571}*cell[12] + x158 + x171 - x172 + x176 + x177 + x180 + x181) - x162*(cell[12] + V{0.0321428571428571});
auto x13 = cell[13] + x161*(-V{0.871428571428571}*cell[13] + x145 - x148 + x166 + x181 + x182 + x183) - x162*(cell[13] + V{0.0321428571428571});
auto x14 = cell[14] + x43*(-V{0.904761904761905}*cell[14] + x185 + x186 + x42 + x48 + x55) - x44*(cell[14] + V{0.0476190476190476});
auto x15 = cell[15] + x43*(-V{0.904761904761905}*cell[15] + x185 + x187 + x38 + x41 + x52 + x53 + x54) - x44*(cell[15] + V{0.0476190476190476});
auto x16 = cell[16] + x43*(-V{0.904761904761905}*cell[16] + x184 + x186 + x187 + x29 + x34 + x40 + x45 + x46 + x56) - x44*(cell[16] + V{0.0476190476190476});
auto x17 = cell[17] + x82*(-V{0.885714285714286}*cell[17] + x115 + x131 + x188 + x190 + x192 + x81 + x87 - x98) - x83*(cell[17] + V{0.0380952380952381});
auto x18 = cell[18] + x82*(-V{0.885714285714286}*cell[18] + x126 + x130 + x137 + x190 + x61 + x65 - x80 + x96 + x99) - x83*(cell[18] + V{0.0380952380952381});
auto x19 = cell[19] + x82*(-V{0.885714285714286}*cell[19] + x108 + x113 - x117 + x136 + x194 + x195 + x62 + x73 + x90) - x83*(cell[19] + V{0.0380952380952381});
auto x20 = cell[20] + x82*(-V{0.885714285714286}*cell[20] - x107 + x118 + x125 + x134 + x191 + x195 + x196 + x66 + x86) - x83*(cell[20] + V{0.0380952380952381});
auto x21 = cell[21] + x82*(-V{0.885714285714286}*cell[21] + x109 + x129 + x130 + x132 - x133 + x194 + x197 + x91 + x93 + x94) - x83*(cell[21] + V{0.0380952380952381});
auto x22 = cell[22] + x82*(-V{0.885714285714286}*cell[22] + x111 - x128 + x135 + x192 + x196 + x197 + x70 + x88 + x89) - x83*(cell[22] + V{0.0380952380952381});
auto x23 = cell[23] + x161*(-V{0.871428571428571}*cell[23] - x167 + x168 + x173 + x177 + x183 + x199 + x200 + x201) - x162*(cell[23] + V{0.0321428571428571});
auto x24 = cell[24] + x161*(-V{0.871428571428571}*cell[24] + x138 - x142 + x143 + x150 + x170 + x176 + x199 + x202) - x162*(cell[24] + V{0.0321428571428571});
auto x25 = cell[25] + x161*(-V{0.871428571428571}*cell[25] + x142 - x143 + x148 + x154 + x170 + x174 + x182 + x200 + x203) - x162*(cell[25] + V{0.0321428571428571});
auto x26 = cell[26] + x161*(-V{0.871428571428571}*cell[26] + x153 + x160 + x163 + x167 - x168 + x172 + x180 + x201 + x202 + x203) - x162*(cell[26] + V{0.0321428571428571});
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[19] + cell[1] + cell[20] + cell[21] + cell[22] + cell[23] + cell[24] + cell[25] + cell[26] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1}, V{0} };
}
};

}

}
