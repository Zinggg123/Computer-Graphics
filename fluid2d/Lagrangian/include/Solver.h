#pragma once
#ifndef __LAGRANGIAN_2D_SOLVER_H__
#define __LAGRANGIAN_2D_SOLVER_H__

#include "ParticleSystem2d.h"
#include "Configure.h"

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        class Solver
        {
        public:
            Solver(ParticleSystem2d &ps);

            void solve();

            void calacc(ParticleInfo2d& p);
            void calx(ParticleInfo2d& p);

        private:
            ParticleSystem2d &mPs;
        };
    }
}

#endif
