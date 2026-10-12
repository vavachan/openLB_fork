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
struct CSE<BoolakeeLinearElasticity<T, descriptors::D2Q8<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x37 = parameters.template get<descriptors::OMEGA_SOLID>()[5];
auto x32 = parameters.template get<descriptors::OMEGA_SOLID>()[0];
auto x34 = parameters.template get<descriptors::OMEGA_SOLID>()[2];
auto x35 = parameters.template get<descriptors::OMEGA_SOLID>()[3];
auto x33 = parameters.template get<descriptors::OMEGA_SOLID>()[1];
auto x18 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x19 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x36 = parameters.template get<descriptors::OMEGA_SOLID>()[4];
auto x8 = V{0.5}*x18;
auto x9 = V{1}*cell[0];
auto x10 = V{0.666666666666667}*cell[5];
auto x11 = V{0.666666666666667}*cell[7];
auto x12 = -V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[6];
auto x13 = x35*(V{0.333333333333333}*cell[0] - V{0.333333333333333}*cell[2] + x10 - x11 + x12 + V{0.166666666666667}*x18);
auto x14 = V{0.5}*x13;
auto x15 = -cell[3];
auto x16 = x33*(cell[0] - cell[1] + cell[2] + x15);
auto x17 = V{0.25}*x16;
auto x20 = -x17;
auto x21 = cell[4] + cell[6];
auto x22 = cell[5] + cell[7];
auto x23 = x21 + x22;
auto x24 = x23*x37;
auto x25 = V{0.5}*x24;
auto x26 = V{1}*cell[1];
auto x27 = V{1}*cell[3];
auto x28 = x26 + x27;
auto x29 = V{1}*cell[2];
auto x30 = x29 + x9;
auto x31 = V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x28 + x30;
auto x38 = x31*x34;
auto x39 = -V{0.25}*x38;
auto x40 = x25 + x39;
auto x41 = x20 + x40;
auto x42 = V{0.5}*x19;
auto x43 = x36*(V{0.333333333333333}*cell[1] - V{0.333333333333333}*cell[3] - x10 + x11 + x12 + V{0.166666666666667}*x19);
auto x44 = V{0.5}*x43;
auto x45 = x17 + x40;
auto x46 = V{1}*cell[4];
auto x47 = V{0.25}*x13;
auto x48 = -cell[5];
auto x49 = -cell[7];
auto x50 = x21 + x48 + x49;
auto x51 = x32*x50;
auto x52 = V{0.25}*x51;
auto x53 = V{0.25}*x43;
auto x54 = V{0.25}*x24;
auto x55 = -x54;
auto x56 = x53 + x55;
auto x57 = V{1}*cell[5];
auto x58 = V{1}*cell[6];
auto x59 = -x58;
auto x60 = x47 + x52;
auto x61 = V{1}*cell[7];
auto x62 = -x29;
auto x63 = x59 + x61;
auto x64 = x46 - x57 + x63;
auto x65 = V{0.75}*x18 + x62 + x64 + x9;
auto x66 = x46 + x57;
auto x67 = x59 - x61 + x66;
auto x68 = V{0.75}*x19 + x26 - x27 + x67;
auto x69 = x58 + x61 + x66;
auto x70 = -cell[2];
auto x71 = -cell[6];
auto x72 = cell[4] + x71;
auto x73 = x39 + x69;
auto x74 = V{1}*x50*(V{0.5}*x32 + V{-1});
cell[0] = -x14 + x41 + x8 + x9;
cell[1] = x26 + x42 - x44 + x45;
cell[2] = x14 + x29 + x41 - x8;
cell[3] = x27 - x42 + x44 + x45;
cell[4] = x46 + x47 - x52 + x56;
cell[5] = -x47 + x52 + x56 + x57;
cell[6] = -x53 - x54 - x59 - x60;
cell[7] = -x53 + x55 + x60 + x61;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[0] = x65;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[1] = x68;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[2] = x46 - V{0.5}*x51 - x57 - x63;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[3] = x31 - V{0.5}*x38;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[4] = -V{0.5}*x16 - x28 - x62 + x9;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[5] = x14 + x64;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[6] = x44 + x67;
cell.template getFieldPointer<olb::descriptors::BARED_MOMENT_VECTOR>()[7] = -x25 + x69;
cell.template getFieldPointer<olb::descriptors::DISP_SOLID>()[0] = x65;
cell.template getFieldPointer<olb::descriptors::DISP_SOLID>()[1] = x68;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[0] = V{1}*cell[0] + V{1}*cell[7] + V{1}*x18 + V{1}*x48 + V{1}*x70 + V{1}*x72;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[1] = V{1}*cell[1] + V{1}*cell[5] + V{1}*x15 + V{1}*x19 + V{1}*x49 + V{1}*x72;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[2] = V{1}*cell[4] - V{1}*x22 - V{1}*x51 - V{1}*x71;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[3] = x31 - x38;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[4] = V{1}*cell[0] - V{1}*cell[1] - V{1}*cell[3] - V{1}*x16 - V{1}*x70;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[5] = x13 + x64;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[6] = x43 + x67;
cell.template getFieldPointer<olb::descriptors::MOMENT_VECTOR>()[7] = V{1}*x23 - V{1}*x24;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[0] = -x20 - x30 - x73;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[1] = x74;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[2] = x74;
cell.template getFieldPointer<olb::descriptors::SIGMA_SOLID>()[3] = -x17 - x28 - x73;
return { V{-1}, V{-1} };
}
};

}

}
