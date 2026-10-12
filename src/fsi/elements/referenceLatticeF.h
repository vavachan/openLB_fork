/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef FSI_REFERENCE_LATTICE_F_H
#define FSI_REFERENCE_LATTICE_F_H

#include "fsi/fields.h"

#include "utilities/geometricOperations.h"

namespace olb {

struct ReferenceLatticePorosityF {
  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_LOWER,
    fields::fsi::ELEMENT_UPPER,
    fields::fsi::ELEMENT_ROTATION,
    fields::fsi::ELEMENT_ROTATION_MATRIX,
    fields::fsi::ELEMENT_U_TRANSLATION,
    fields::fsi::ELEMENT_U_ROTATION,
    fields::fsi::ELEMENT_REFERENCE_DELTA_X,
    fields::fsi::ELEMENT_REFERENCE_PROJECTION,
    fields::fsi::ELEMENT_REFERENCE_POROSITY,
    fields::fsi::PERIODIC_DIRECTIONS,
    fields::fsi::DOMAIN_LENGTH
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_LOWER>,
    fields::array_of<fields::fsi::ELEMENT_UPPER>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_ROTATION>,
    fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    fields::array_of<fields::fsi::ELEMENT_U_ROTATION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITY>,
    fields::array_of<fields::fsi::PERIODIC_DIRECTIONS>,
    fields::array_of<fields::fsi::DOMAIN_LENGTH>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <typename PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned iElement) any_platform {
    return params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>()[iElement];
  }

  template <typename PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto lowerBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_LOWER>>();
    const auto upperBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_UPPER>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION>>();
    const auto rotationMatricies = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>>();

    const auto periodicDirsArr = params.template get<fields::array_of<fields::fsi::PERIODIC_DIRECTIONS>>();
    const auto domainLengthsArr = params.template get<fields::array_of<fields::fsi::DOMAIN_LENGTH>>();
    const Vector<V, DESCRIPTOR::d> periodicDirs{periodicDirsArr,iElement};
    const Vector<V, DESCRIPTOR::d> domainLengths{domainLengthsArr,iElement};

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION> elementRotation{rotations, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION_MATRIX> elementRotationMatrix{rotationMatricies, iElement};


    //std::cout<<periodicDirs<<"\t"<<domainLengths<<"\n";
    Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    
    for(int d=0; d<DESCRIPTOR::d; d++)
    {
	    if(periodicDirs[d])
		    shiftedPhysR[d] -= domainLengths[d]*std::round(shiftedPhysR[d]/domainLengths[d]);
    }

    //const auto rotatedPhysR = util::matrixVectorProduct(util::calculateRotationMatrix<V,DESCRIPTOR::d>(-1*elementRotation),
                                                        //shiftedPhysR);
    auto inverseRotationMatrix = util::invertRotationMatrix<V,DESCRIPTOR::d>(elementRotationMatrix);

    const auto rotatedPhysR = util::matrixVectorProduct(inverseRotationMatrix,shiftedPhysR);
    const Vector<V,DESCRIPTOR::d> shiftedLowerR = elementLowerR - elementPivotR;
    const Vector<V,DESCRIPTOR::d> shiftedUpperR = elementUpperR - elementPivotR;

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    //std::cout<<elementPivotR<<"\n";
    //std::cout<<shiftedLowerR<<"\t"<<shiftedUpperR<<"\t"<<rotatedPhysR<<"\t"<<physDeltaX<<"\n";
    //std::cout<<elementLowerR<<"\t"<<elementUpperR<<"\t"<<elementPivotR<<"\t"<<physDeltaX<<"\n";
    if (shiftedLowerR+physDeltaX < rotatedPhysR && rotatedPhysR < shiftedUpperR-physDeltaX) {
      const auto refDeltaX = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>();
      const auto refProjection = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>();
      const auto refLattice = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITY>>();

      const Vector<std::size_t,DESCRIPTOR::d> projection{refProjection, iElement};
      const Vector<int,DESCRIPTOR::d> discreteR = round((rotatedPhysR - shiftedLowerR) / refDeltaX[iElement]);


      //std::cout<<projection<<std::endl;
      auto porosity = refLattice[iElement][projection * discreteR];
      //std::cout<<rotatedPhysR[0]<<"\t"<<rotatedPhysR[1]<<"\t"<<rotatedPhysR[2]<<"\n";
      //std::cout<<physR[0]<<"\t"<<physR[1]<<"\t"<<physR[2]<<"\n";
    //if(!porosity)
    //{
    //        //std::cout<<rotatedPhysR[0]<<"\t"<<rotatedPhysR[1]<<"\t"<<rotatedPhysR[2]<<"\n";
    //        std::cout<<physR[0]<<"\t"<<physR[1]<<"\t"<<physR[2]<<"\n";
    //}

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        porosity += refLattice[iElement][projection * (discreteR + descriptors::c<DESCRIPTOR>(iPop))];
      }
      porosity /= DESCRIPTOR::q;
      return porosity;
    } else {
      return V{1};
    }
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto translationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>>();
    const auto rotationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_ROTATION>>();
    const auto periodicDirsArr = params.template get<fields::array_of<fields::fsi::PERIODIC_DIRECTIONS>>();
    const auto domainLengthsArr = params.template get<fields::array_of<fields::fsi::DOMAIN_LENGTH>>();

    const Vector<V, DESCRIPTOR::d> periodicDirs{periodicDirsArr,iElement};
    const Vector<V, DESCRIPTOR::d> domainLengths{domainLengthsArr,iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};

    Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;


    for(int d=0; d<DESCRIPTOR::d; d++)
    {
	    if(periodicDirs[d])
		    shiftedPhysR[d] -= domainLengths[d]*std::round(shiftedPhysR[d]/domainLengths[d]);
    }

    Vector<V,DESCRIPTOR::d> u{translationUs, iElement};
    if constexpr (DESCRIPTOR::d == 2) {
      u[0] -= rotationUs[iElement] * shiftedPhysR[1];
      u[1] += rotationUs[iElement] * shiftedPhysR[0];
    } else {
      u += crossProduct3D(Vector<V,DESCRIPTOR::d>{rotationUs, iElement},
                          shiftedPhysR);
    }
    u *= params.template get<fields::converter::PHYS_VELOCITY>();
    return u;
  }
};

struct ReferenceLatticeWithWallModelPorosityF : public ReferenceLatticePorosityF {
  using data = typename ReferenceLatticePorosityF::data::template include<
    fields::fsi::ELEMENT_REFERENCE_Y1
  >;

  using parameters = typename ReferenceLatticePorosityF::parameters::template include<
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_Y1>
  >;

  template <typename PARAMETERS, typename PHYS_R>
  auto computeY1(PARAMETERS& params,
                 PHYS_R& physR,
                 unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto lowerBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_LOWER>>();
    const auto upperBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_UPPER>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION>>();

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION> elementRotation{rotations, iElement};

    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    const auto rotatedPhysR = util::matrixVectorProduct(util::calculateRotationMatrix<V,DESCRIPTOR::d>(-1*elementRotation),
                                                        shiftedPhysR);
    const Vector<V,DESCRIPTOR::d> shiftedLowerR = elementLowerR - elementPivotR;
    const Vector<V,DESCRIPTOR::d> shiftedUpperR = elementUpperR - elementPivotR;

    const auto inverseRotation = util::calculateRotationMatrix<V,DESCRIPTOR::d>(elementRotation);

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    if (shiftedLowerR+physDeltaX < rotatedPhysR && rotatedPhysR < shiftedUpperR-physDeltaX) {
      const auto refDeltaX = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>();
      const auto refProjection = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>();
      const auto refLattice = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_Y1>>();

      const Vector<std::size_t,DESCRIPTOR::d> projection{refProjection, iElement};
      const Vector<int,DESCRIPTOR::d> discreteR = round((rotatedPhysR - shiftedLowerR) / refDeltaX[iElement]);

      const Vector<V*,DESCRIPTOR::d> y1s{refLattice, iElement};

      Vector<V,DESCRIPTOR::d> y1{y1s, projection * discreteR};
      for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
        y1 += Vector<V,DESCRIPTOR::d>{y1s, projection * (discreteR + descriptors::c<DESCRIPTOR>(iPop))};
      }
      y1 /= descriptors::q<DESCRIPTOR>();

      return util::matrixVectorProduct(inverseRotation, y1);
    } else {
      return Vector<V,DESCRIPTOR::d>{};
    }
  }
};

struct ArbitrarilyRotateableReferenceLatticeWithWallModelPorosityF {
  /// Data fields in element store
  using data = meta::list<
    fields::fsi::ELEMENT_LOWER,
    fields::fsi::ELEMENT_UPPER,
    fields::fsi::ELEMENT_ROTATION_MATRIX,
    fields::fsi::ELEMENT_U_TRANSLATION,
    fields::fsi::ELEMENT_U_ROTATION,
    fields::fsi::ELEMENT_REFERENCE_DELTA_X,
    fields::fsi::ELEMENT_REFERENCE_PROJECTION,
    fields::fsi::ELEMENT_REFERENCE_POROSITY,
    fields::fsi::ELEMENT_REFERENCE_Y1
  >;

  /// Parameter fields required for computation
  using parameters = meta::list<
    fields::array_of<fields::fsi::ELEMENT_LOWER>,
    fields::array_of<fields::fsi::ELEMENT_UPPER>,
    fields::array_of<fields::fsi::ELEMENT_PIVOT>,
    fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>,
    fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>,
    fields::array_of<fields::fsi::ELEMENT_U_ROTATION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITY>,
    fields::array_of<fields::fsi::ELEMENT_REFERENCE_Y1>,
    // General non-FSI specific parameters
    fields::converter::PHYS_VELOCITY,
    fields::converter::PHYS_DELTA_X
  >;

  template <typename PARAMETERS, typename PHYS_R>
  int tag(PARAMETERS& params,
          PHYS_R& physR,
          unsigned iElement) any_platform {
    return params.template get<fields::array_of<fields::fsi::ELEMENT_TAG>>()[iElement];
  }

  template <typename PARAMETERS, typename PHYS_R>
  bool isInterior(PARAMETERS& params,
                  PHYS_R& physR,
                  unsigned iElement) any_platform {
    return false;
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto compute(PARAMETERS& params,
               PHYS_R& physR,
               unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto lowerBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_LOWER>>();
    const auto upperBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_UPPER>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>>();

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION_MATRIX> elementRotation{rotations, iElement};

    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    const auto rotatedPhysR = util::matrixVectorProduct(util::inverseRotation(elementRotation),
                                                        shiftedPhysR);
    const Vector<V,DESCRIPTOR::d> shiftedLowerR = elementLowerR - elementPivotR;
    const Vector<V,DESCRIPTOR::d> shiftedUpperR = elementUpperR - elementPivotR;

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    if (shiftedLowerR+physDeltaX < rotatedPhysR && rotatedPhysR < shiftedUpperR-physDeltaX) {
      const auto refDeltaX = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>();
      const auto refProjection = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>();
      const auto refLattice = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_POROSITY>>();

      const Vector<std::size_t,DESCRIPTOR::d> projection{refProjection, iElement};
      const Vector<int,DESCRIPTOR::d> discreteR = round((rotatedPhysR - shiftedLowerR) / refDeltaX[iElement]);

      auto porosity = refLattice[iElement][projection * discreteR];
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        porosity += refLattice[iElement][projection * (discreteR + descriptors::c<DESCRIPTOR>(iPop))];
      }
      porosity /= DESCRIPTOR::q;
      return porosity;
    } else {
      return V{1};
    }
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeU(PARAMETERS& params,
                PHYS_R& physR,
                unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>>();

    const auto translationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_TRANSLATION>>();
    const auto rotationUs = params.template get<fields::array_of<fields::fsi::ELEMENT_U_ROTATION>>();

    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION_MATRIX> elementRotation{rotations, iElement};

    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    const auto rotatedPhysR = util::matrixVectorProduct(util::inverseRotation(elementRotation),
                                                        shiftedPhysR);

    Vector<V,DESCRIPTOR::d> u{translationUs, iElement};
    if constexpr (DESCRIPTOR::d == 2) {
      u[0] -= rotationUs[iElement] * rotatedPhysR[1];
      u[1] += rotationUs[iElement] * rotatedPhysR[0];
    } else {
      u += crossProduct3D(Vector<V,DESCRIPTOR::d>{rotationUs, iElement},
                          rotatedPhysR);
    }

    return util::matrixVectorProduct(elementRotation, u)
         * params.template get<fields::converter::PHYS_VELOCITY>();
  }

  template <typename PARAMETERS, typename PHYS_R>
  auto computeY1(PARAMETERS& params,
                 PHYS_R& physR,
                 unsigned iElement) any_platform {
    using V = typename PARAMETERS::value_t;
    using DESCRIPTOR = typename PARAMETERS::descriptor_t;

    const auto lowerBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_LOWER>>();
    const auto upperBounds = params.template get<fields::array_of<fields::fsi::ELEMENT_UPPER>>();
    const auto pivots = params.template get<fields::array_of<fields::fsi::ELEMENT_PIVOT>>();
    const auto rotations = params.template get<fields::array_of<fields::fsi::ELEMENT_ROTATION_MATRIX>>();

    const Vector<V,DESCRIPTOR::d> elementLowerR{lowerBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementUpperR{upperBounds, iElement};
    const Vector<V,DESCRIPTOR::d> elementPivotR{pivots, iElement};
    const FieldD<V,DESCRIPTOR,fields::fsi::ELEMENT_ROTATION_MATRIX> elementRotation{rotations, iElement};

    const Vector<V,DESCRIPTOR::d> shiftedPhysR = physR - elementPivotR;
    const auto rotatedPhysR = util::matrixVectorProduct(util::inverseRotation(elementRotation),
                                                        shiftedPhysR);
    const Vector<V,DESCRIPTOR::d> shiftedLowerR = elementLowerR - elementPivotR;
    const Vector<V,DESCRIPTOR::d> shiftedUpperR = elementUpperR - elementPivotR;

    const auto inverseRotation = elementRotation;

    const auto physDeltaX = params.template get<fields::converter::PHYS_DELTA_X>();

    if (shiftedLowerR+physDeltaX < rotatedPhysR && rotatedPhysR < shiftedUpperR-physDeltaX) {
      const auto refDeltaX = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_DELTA_X>>();
      const auto refProjection = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_PROJECTION>>();
      const auto refLattice = params.template get<fields::array_of<fields::fsi::ELEMENT_REFERENCE_Y1>>();

      const Vector<std::size_t,DESCRIPTOR::d> projection{refProjection, iElement};
      const Vector<int,DESCRIPTOR::d> discreteR = round((rotatedPhysR - shiftedLowerR) / refDeltaX[iElement]);

      const Vector<V*,DESCRIPTOR::d> y1s{refLattice, iElement};

      Vector<V,DESCRIPTOR::d> y1{y1s, projection * discreteR};
      for (int iPop=1; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
        y1 += Vector<V,DESCRIPTOR::d>{y1s, projection * (discreteR + descriptors::c<DESCRIPTOR>(iPop))};
      }
      y1 /= descriptors::q<DESCRIPTOR>();

      return util::matrixVectorProduct(inverseRotation, y1);
    } else {
      return Vector<V,DESCRIPTOR::d>{};
    }
  }

};

}

#endif
