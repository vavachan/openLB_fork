/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef DERIVATIVE_PRIMITIVE_F_H
#define DERIVATIVE_PRIMITIVE_F_H

#include "../concept.h"

namespace olb {

namespace functors {

struct DissipationDF {
  using parameters = meta::list<descriptors::OMEGA,
                                descriptors::DT,
                                descriptors::PHYS_VISCOSITY>;

  using result_t = descriptors::FIELD_MATRIX<descriptors::DISSIPATION,
                                         descriptors::POPULATION>;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V dt = parameters.template get<descriptors::DT>();
    const V physViscosity = parameters.template get<descriptors::PHYS_VISCOSITY>();

    V rho = 0; V u[DESCRIPTOR::d]{0}; V pi[util::TensorVal<DESCRIPTOR>::n]{0};
    cell.computeAllMomenta(rho, u, pi);

    FieldD<V,DESCRIPTOR,result_t> dissipationD{};
    for (std::size_t jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      dissipationD[jPop] = 0;
      V dpidf[util::TensorVal<DESCRIPTOR>::n];
      opti::dualLbMomentaHelpers<DESCRIPTOR>::dPiDf(cell, dpidf, jPop);
      for (std::size_t iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
        for (std::size_t iBeta=0; iBeta < DESCRIPTOR::d; ++iBeta) {
          const int iPi = util::serialSymmetricTensorIndex<DESCRIPTOR::d>(iAlpha, iBeta);
          dissipationD[jPop] += pi[iPi] * dpidf[iPi] - pi[iPi] * pi[iPi] / rho;
        }
      }
      dissipationD[jPop] *= 2. * util::pow(omega*descriptors::invCs2<V,DESCRIPTOR>() / rho, 2)
                            / 2. * physViscosity / dt / dt;
    }
    return dissipationD;
  }
};

struct PorousDissipationDF {
  using parameters = meta::list<descriptors::OMEGA,
                                descriptors::DX,
                                descriptors::VISCOSITY,
                                descriptors::CONVERSION_VELOCITY,
                                descriptors::PHYS_VISCOSITY>;

  using result_t = descriptors::FIELD_MATRIX<descriptors::POROUS_DISSIPATION,
                                         descriptors::POPULATION>;
  using fields_t = meta::list<descriptors::POROSITY>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V dx = parameters.template get<descriptors::DX>();
    const V viscosity = parameters.template get<descriptors::VISCOSITY>();
    const V physViscosity = parameters.template get<descriptors::PHYS_VISCOSITY>();
    const V conversionVelocity = parameters.template get<descriptors::CONVERSION_VELOCITY>();

    Vector<V,DESCRIPTOR::d> u{0};
    cell.computeU(u.data());
    const V gridTerm = dx * dx * viscosity / omega;
    const V invPermeability = (V(1) - porosity) / gridTerm;
    Vector<V,DESCRIPTOR::d> dudf;
    FieldD<V,DESCRIPTOR,result_t> dissipationD{};
    for (std::size_t jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
      opti::dualLbMomentaHelpers<DESCRIPTOR>::dUDf(cell, dudf, jPop);
      dissipationD[jPop] = physViscosity * invPermeability * 2. * (u * dudf) *
                       conversionVelocity * conversionVelocity;
    }
    return dissipationD;
  }
};

struct TotalDissipationDalpha {
  using parameters = meta::list<descriptors::OMEGA,
                                descriptors::DX,
                                descriptors::VISCOSITY,
                                descriptors::CONVERSION_VELOCITY,
                                descriptors::PHYS_VISCOSITY,
                opti::REG_ALPHA>;

  using result_t = descriptors::FIELD_MATRIX<descriptors::DISSIPATION,
                                         descriptors::POROSITY>;
  using fields_t = meta::list<descriptors::POROSITY>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V porosity = cell.template getField<descriptors::POROSITY>();
    const V omega = parameters.template get<descriptors::OMEGA>();
    const V dx = parameters.template get<descriptors::DX>();
    const V viscosity = parameters.template get<descriptors::VISCOSITY>();
    const V physViscosity = parameters.template get<descriptors::PHYS_VISCOSITY>();
    const V conversionVelocity = parameters.template get<descriptors::CONVERSION_VELOCITY>();
    const V regAlpha = parameters.template get<opti::REG_ALPHA>();

    Vector<V,DESCRIPTOR::d> u{0};
    cell.computeU(u.data());
    const V gridTerm = dx * dx * viscosity / omega;
    const V invPermeability = (V(1) - porosity) / gridTerm;
    const V uNormSq = util::euklidN2(u.data(), DESCRIPTOR::d) * conversionVelocity
                  * conversionVelocity;
    FieldD<V,DESCRIPTOR,result_t> dJDalpha = physViscosity * uNormSq * invPermeability + regAlpha;
    return dJDalpha;
  }
};

}

}

#endif
