#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"

namespace FluidSimulation
{
    namespace Eulerian3d
    {
        Solver::Solver(MACGrid3d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::solve()
        {
            add_forces();
            advect();
            project();
            extrapolate();
            advect_dens();
        }

        void Solver::add_forces() {
            FOR_EACH_CELL
            {
                double bou = mGrid.getBoussinesqForce(glm::vec3((i + 0.5) * mGrid.cellSize, (j + 0.5) * mGrid.cellSize, (k + 0.5) * mGrid.cellSize));
                mGrid.mW(i, j, k) += (bou)*Eulerian3dPara::dt;
            }
        }

        void Solver::advect() {
            Glb::GridData3dX mU_prev(mGrid.mU);
            Glb::GridData3dY mV_prev(mGrid.mV);
            Glb::GridData3dZ mW_prev(mGrid.mW);
            float dt0 = Eulerian3dPara::dt;
            FOR_EACH_FACE
            {
                if (!mGrid.isSolidCell(i - 1, j, k) && j < mGrid.dim[1] && k < mGrid.dim[2]) {
                    glm::vec3 tmp1((i)*mGrid.cellSize, (j + 0.5) * mGrid.cellSize, (k + 0.5) * mGrid.cellSize);
                    glm::vec3 p1 = mGrid.semiLagrangian(tmp1, dt0);
                    mU_prev(i, j, k) = mGrid.mU.interpolate(p1);
                }

                if (!mGrid.isSolidCell(i, j - 1, k) && i < mGrid.dim[0] && k < mGrid.dim[2]) {
                    glm::vec3 tmp2((i + 0.5) * mGrid.cellSize, (j)*mGrid.cellSize, (k + 0.5) * mGrid.cellSize);
                    glm::vec3 p2 = mGrid.semiLagrangian(tmp2, dt0);
                    mV_prev(i, j, k) = mGrid.mV.interpolate(p2);
                }

                if (!mGrid.isSolidCell(i, j, k - 1) && i < mGrid.dim[0] && j < mGrid.dim[1]) {
                    glm::vec3 tmp3((i + 0.5) * mGrid.cellSize, (j + 0.5) * mGrid.cellSize, (k)*mGrid.cellSize);
                    glm::vec3 p3 = mGrid.semiLagrangian(tmp3, dt0);
                    mW_prev(i, j, k) = mGrid.mW.interpolate(p3);
                }
            }
            mGrid.mU = mU_prev;
            mGrid.mV = mV_prev;
            mGrid.mW = mW_prev;

        }

        void Solver::project() {
            int nx = mGrid.dim[0]; 
            int ny = mGrid.dim[1]; 
            int nz = mGrid.dim[2]; 

            std::vector<std::vector<std::vector<double>>> p(nx + 2,
                std::vector<std::vector<double>>(ny + 2,
                    std::vector<double>(nz + 2, 0.0)));
            std::vector<std::vector<std::vector<double>>> div(nx,
                std::vector<std::vector<double>>(ny,
                    std::vector<double>(nz, 0.0)));

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        div[i][j][k] = mGrid.getDivergence(i, j, k);
                        if (j == 0 || i == 0 || k == 0) { 
                            div[i][j][k] = (mGrid.mU(i + 1, j, k) - mGrid.mU(i, j, k) + mGrid.mV(i, j + 1, k) - mGrid.mV(i, j, k)
                                + mGrid.mW(i, j, k + 1) - mGrid.mW(i, j, k)) / mGrid.cellSize;
                        }
                    }
                }
            }
            for (int iter = 0; iter < 200; ++iter) {
                for (int i = 1; i <= nx; ++i) {
                    for (int j = 1; j <= ny; ++j) {
                        for (int k = 1; k <= nz; ++k) {
                            p[i][j][k] = 1.0 / 6.0 * (p[i + 1][j][k] + p[i - 1][j][k] +
                                p[i][j + 1][k] + p[i][j - 1][k] +
                                p[i][j][k + 1] + p[i][j][k - 1] -
                                div[i - 1][j - 1][k - 1]);
                        }
                    }
                }
                for (int i = 1; i <= nx; ++i) {
                    for (int j = 1; j <= ny; ++j) {
                        p[i][j][0] = p[i][j][1];
                        p[i][j][nz + 1] = p[i][j][nz];
                    }
                }
                for (int i = 1; i <= nx; ++i) {
                    for (int k = 1; k <= nz; ++k) {
                        p[i][0][k] = p[i][1][k];
                        p[i][ny + 1][k] = p[i][ny][k];
                    }
                }
                for (int j = 1; j <= ny; ++j) {
                    for (int k = 1; k <= nz; ++k) {
                        p[0][j][k] = p[1][j][k];
                        p[nx + 1][j][k] = p[nx][j][k];
                    }
                }
            }
            for (int i = 1; i <= nx; ++i) {
                for (int j = 1; j <= ny; ++j) {
                    for (int k = 1; k <= nz; ++k) {
                        if (!mGrid.isSolidCell(i - 1, j - 1, k - 1)) {
                            mGrid.mU(i - 1, j - 1, k - 1) -= 0.5 * (p[i][j][k] - p[i - 1][j][k]) / mGrid.cellSize;
                            mGrid.mV(i - 1, j - 1, k - 1) -= 0.5 * (p[i][j][k] - p[i][j - 1][k]) / mGrid.cellSize;
                            mGrid.mW(i - 1, j - 1, k - 1) -= 0.5 * (p[i][j][k] - p[i][j][k - 1]) / mGrid.cellSize;
                        }
                    }
                }
            }
        }

        void Solver::extrapolate() {
            int nx = mGrid.dim[0]; 
            int ny = mGrid.dim[1]; 
            int nz = mGrid.dim[2]; 

            for (int j = 0; j <= ny; ++j) {
                for (int k = 0; k <= nz; ++k) {
                    mGrid.mU(0, j, k) = 0.0;
                    mGrid.mU(nx, j, k) = 0.0;
                }
            }
            for (int i = 0; i <= nx; ++i) {
                for (int k = 0; k <= nz; ++k) {
                    mGrid.mV(i, 0, k) = 0.0;
                    mGrid.mV(i, ny, k) = 0.0;
                }
            }
            for (int i = 0; i <= nx; ++i) {
                for (int j = 0; j <= ny; ++j) {
                    mGrid.mW(i, j, 0) = 0.0;
                    mGrid.mW(i, j, nz) = 0.0;
                }
            }
        }

        void Solver::advect_dens() {
            Glb::CubicGridData3d mD_prev(mGrid.mD);
            Glb::CubicGridData3d mT_prev(mGrid.mT);

            float dt0 = Eulerian3dPara::dt;
            FOR_EACH_CELL
            {
                glm::vec3 tmp((i + 0.5) * mGrid.cellSize, (j + 0.5) * mGrid.cellSize, (k + 0.5) * mGrid.cellSize);
                glm::vec3 p = mGrid.semiLagrangian(tmp, dt0);
                mD_prev(i, j, k) = mGrid.mD.interpolate(p);
                mT_prev(i, j, k) = mGrid.mT.interpolate(p);
            }
            mGrid.mD = mD_prev;
            mGrid.mT = mT_prev;

        }
    }
}
