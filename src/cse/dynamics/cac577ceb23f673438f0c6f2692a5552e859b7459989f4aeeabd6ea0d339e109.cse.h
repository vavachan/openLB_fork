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
struct CSE<BoolakeeLinearElasticityBoundary<T, descriptors::D2Q8<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(7);
auto x56 = parameters.template get<descriptors::OMEGA_SOLID>()[0];
auto x58 = parameters.template get<descriptors::OMEGA_SOLID>()[2];
auto x16 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(0);
auto x60 = parameters.template get<descriptors::OMEGA_SOLID>()[4];
auto x57 = parameters.template get<descriptors::OMEGA_SOLID>()[1];
auto x19 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(3);
auto x18 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(2);
auto x17 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(1);
auto x59 = parameters.template get<descriptors::OMEGA_SOLID>()[3];
auto x61 = parameters.template get<descriptors::OMEGA_SOLID>()[5];
auto x21 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(5);
auto x27 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x20 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(4);
auto x22 = cell.template getFieldComponent<olb::descriptors::CURRENT_CELL>(6);
auto x26 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x8 = V{0.666666666666667}*cell[5];
auto x9 = V{0.666666666666667}*cell[7];
auto x10 = -V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[6];
auto x11 = x59*(V{0.333333333333333}*cell[0] - V{0.333333333333333}*cell[2] + x10 + x8 - x9);
auto x12 = V{0.5}*x11;
auto x13 = V{0.25}*cell[2];
auto x14 = V{0.5}*x26;
auto x15 = x58*(cell[0] - cell[1] + cell[2] - cell[3]);
auto x24 = V{0.25}*x15;
auto x25 = V{0.25}*cell[1];
auto x28 = V{0.25}*cell[3];
auto x29 = V{1} - 1/x57;
auto x30 = V{1} / (x29);
auto x31 = V{0.0625}*x30;
auto x32 = V{0.25} - x31;
auto x33 = V{1}*cell[1];
auto x34 = V{1}*cell[3];
auto x35 = x33 + x34;
auto x36 = V{1}*cell[0];
auto x37 = V{1}*cell[2];
auto x38 = x36 + x37;
auto x39 = V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x35 + x38;
auto x40 = x39*x57;
auto x41 = x39 - x40;
auto x42 = V{1} - V{0.25}*x30;
auto x43 = V{0.5}*cell[4]*x42;
auto x44 = V{0.5}*cell[5]*x42;
auto x45 = V{0.5}*cell[6]*x42;
auto x46 = V{0.5}*cell[7]*x42;
auto x47 = -1/x29;
auto x48 = V{0.25}*x47 + V{1};
auto x49 = V{0.125}*x47;
auto x50 = x61*(cell[0]*x49 + cell[1]*x49 + cell[2]*x49 + cell[3]*x49 + cell[4]*x48 + cell[5]*x48 + cell[6]*x48 + cell[7]*x48);
auto x51 = V{0.5}*x50;
auto x52 = -V{0.0625}*cell[0]*x30 - V{0.0625}*cell[1]*x30 - V{0.0625}*cell[2]*x30 - V{0.0625}*cell[3]*x30 + x24 + x25 + x28 - x32*x41 + x43 + x44 + x45 + x46 - x51;
auto x53 = V{0.75}*cell[0] - x12 - x13 + x14 - x52;
auto x54 = V{0.5}*x27;
auto x55 = x60*(V{0.333333333333333}*cell[1] - V{0.333333333333333}*cell[3] + x10 - x8 + x9);
auto x62 = V{0.5}*x55;
auto x63 = V{0.25}*cell[0];
auto x64 = cell[0]*x31 + cell[1]*x31 + cell[2]*x31 + cell[3]*x31 - x13 + x24 + x32*x41 - x43 - x44 - x45 - x46 + x51 - x63;
auto x65 = V{0.75}*cell[1] - x28 + x54 - x62 + x64;
auto x66 = V{0.75}*cell[2] + x12 - x14 - x52 - x63;
auto x67 = V{0.75}*cell[3] - x25 - x54 + x62 + x64;
auto x68 = V{0.25}*x11;
auto x69 = -x68;
auto x70 = V{0.25}*cell[6];
auto x71 = V{0.25}*cell[7];
auto x72 = V{0.03125}*cell[0];
auto x73 = V{0.03125}*cell[1];
auto x74 = V{0.03125}*cell[2];
auto x75 = V{0.03125}*cell[3];
auto x76 = -x30*(V{0.0625}*cell[4] + V{0.0625}*cell[5] + V{0.0625}*cell[6] + V{0.0625}*cell[7] - V{0.03125}*x40 + x72 + x73 + x74 + x75);
auto x77 = x30*x72;
auto x78 = x30*x73;
auto x79 = x30*x74;
auto x80 = x30*x75;
auto x81 = V{0.25}*x50;
auto x82 = -V{0.25}*cell[4]*x42;
auto x83 = -V{0.25}*cell[5]*x42;
auto x84 = -V{0.25}*cell[6]*x42;
auto x85 = -V{0.25}*cell[7]*x42;
auto x86 = V{0.25}*x55;
auto x87 = x70 + x71 + x76 + x77 + x78 + x79 + x80 + x81 + x82 + x83 + x84 + x85 - x86;
auto x88 = V{0.25}*cell[5];
auto x89 = cell[4] - cell[5] + cell[6] - cell[7];
auto x90 = x56*x89;
auto x91 = V{0.25}*x90;
auto x92 = x88 + x91;
auto x93 = V{0.75}*cell[4] - x69 - x87 - x92;
auto x94 = -x91;
auto x95 = V{0.25}*cell[4];
auto x96 = x68 + x95;
auto x97 = V{0.75}*cell[5] - x87 - x94 - x96;
auto x98 = x76 + x77 + x78 + x79 + x80 + x81 + x82 + x83 + x84 + x85 + x86;
auto x99 = V{0.75}*cell[6] - x71 - x92 - x96 - x98;
auto x100 = V{0.75}*cell[7] - x69 - x70 - x88 - x94 - x95 - x98;
auto x101 = -x37;
auto x102 = V{1}*cell[4];
auto x103 = V{1}*cell[5];
auto x104 = V{1}*cell[7];
auto x105 = V{1}*cell[6];
auto x106 = -x105;
auto x107 = x104 + x106;
auto x108 = x102 - x103 + x107;
auto x109 = x101 + x108 + x36;
auto x110 = x109 + x14;
auto x111 = x102 + x103;
auto x112 = -x104 + x106 + x111;
auto x113 = x112 + x33 - x34;
auto x114 = x113 + x54;
auto x115 = V{0.125}*x30;
auto x116 = cell[0]*x115 + cell[1]*x115 + cell[2]*x115 + cell[3]*x115;
auto x117 = x104 + x105 + x111 - V{0.25}*x39*x57;
auto x118 = V{1}*x89*(V{0.5}*x56 + V{-1});
cell[0] = x53;
cell[1] = x65;
cell[2] = x66;
cell[3] = x67;
cell[4] = x93;
cell[5] = x97;
cell[6] = x99;
cell[7] = x100;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[0] = x110;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[1] = x114;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[2] = x102 - x103 - x107 - V{0.5}*x90;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[3] = x39 - V{0.5}*x40;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[4] = -x101 - V{0.5}*x15 - x35 + x36;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[5] = x108 + x12;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[6] = x112 + x62;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[7] = V{1}*cell[4]*x42 + V{1}*cell[5]*x42 + V{1}*cell[6]*x42 + V{1}*cell[7]*x42 - x116 - x51;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[0] = x53;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[1] = x65;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[2] = x66;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[3] = x67;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[4] = x93;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[5] = x97;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[6] = x99;
cell.template getFieldPointer<olb::descriptors::CURRENT_CELL>()[7] = x100;
cell.template getFieldPointer<olb::descriptors::DISP_SOLID>()[0] = x110;
cell.template getFieldPointer<olb::descriptors::DISP_SOLID>()[1] = x114;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[0] = x109 + x26;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[1] = x113 + x27;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[2] = V{1}*cell[4] - V{1}*cell[5] + V{1}*cell[6] - V{1}*cell[7] - V{1}*x90;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[3] = x41;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[4] = V{1}*cell[0] - V{1}*cell[1] + V{1}*cell[2] - V{1}*cell[3] - V{1}*x15;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[5] = x108 + x11;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[6] = x112 + x55;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[7] = cell[4]*x42 + cell[5]*x42 + cell[6]*x42 + cell[7]*x42 - x116 - x50;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[0] = x16;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[1] = x17;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[2] = x18;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[3] = x19;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[4] = x20;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[5] = x21;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[6] = x22;
cell.template getFieldPointer<olb::descriptors::PREVIOUS_CELL>()[7] = x23;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[0] = -x117 + x24 - x38;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[1] = x118;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[2] = x118;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[3] = -x117 - x24 - x35;
return { V{-1}, V{-1} };
}
};

}

}
