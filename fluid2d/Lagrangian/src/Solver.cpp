#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        Solver::Solver(ParticleSystem2d& ps) : mPs(ps)
        {
        }


        //void Solver::solve()
        //{
        //    // TODO
        //    // Solves the fluid simulation by performing some steps, which may include:
        //    // 1. compute density 
        //    // 2. compute press
        //    // 3. compute accleration
        //    // 4. update velocity and position
        //    // 5. check boundary
        //    // 6. update block id
        //    // ...

        //    std::cout << "1ps: " << mPs.particles.size() << std::endl;
        //    calpre();
        //    std::cout << "2ps: " << mPs.particles.size() << std::endl;
        //    calacc();
        //    calx();
        //    mPs.updateBlockInfo();
        //}


    // SolverÀàµÄsolveº¯Êý
        void Solver::solve() {
            float deltaTime = Lagrangian2dPara::dt;
            int i = 0;
            for (auto& p : mPs.particles) {
                calacc(p); //std::cout << i <<" ";
                calx(p);

                p.accleration += glm::vec2(Lagrangian2dPara::gravityX, -Lagrangian2dPara::gravityY);
                p.velocity += p.accleration * deltaTime;
                p.position += p.velocity * deltaTime;

                if (p.position.x < mPs.lowerBound.x + Lagrangian2dPara::eps) {
                    p.position.x = mPs.lowerBound.x + Lagrangian2dPara::eps;
                    p.velocity.x *= -Lagrangian2dPara::velocityAttenuation;
                }
                if (p.position.y < mPs.lowerBound.y + Lagrangian2dPara::eps) {
                    p.position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                    p.velocity.y *= -Lagrangian2dPara::velocityAttenuation;
                }
                if (p.position.x > mPs.upperBound.x + Lagrangian2dPara::eps) {
                    p.position.x = mPs.upperBound.x + Lagrangian2dPara::eps;
                    p.velocity.x *= -Lagrangian2dPara::velocityAttenuation;
                }
                if (p.position.y > mPs.upperBound.y + Lagrangian2dPara::eps) {
                    p.position.y = mPs.upperBound.y + Lagrangian2dPara::eps;
                    p.velocity.y *= -Lagrangian2dPara::velocityAttenuation;
                }

                if (length(p.velocity) > Lagrangian2dPara::maxVelocity) {
                    p.velocity = normalize(p.velocity) * Lagrangian2dPara::maxVelocity;
                }
                i++;
            }
            mPs.updateBlockInfo();
        }

        void Solver::calacc(ParticleInfo2d& p) {
            int bid;
            glm::uvec2 ran;
            float rhosum = 0; 
            float mploy = 315.0 / (64.0 * 3.14159f * pow(mPs.supportRadius, 9));
            int sz = int(sqrt(mPs.blockExtens.size()));
            for (int off : mPs.blockIdOffs) {
                bid = p.blockId + off;
                if (p.blockId % sz == 0 && (bid + 1) % sz == 0 || (p.blockId+1) % sz == 0 && bid % sz == 0) {
                    continue;
                }
                ran = mPs.blockExtens[bid];

                if (ran.x == ran.y)continue;
                for (int i = ran.x; i < ran.y; i++) {
                    ParticleInfo2d& p2 = mPs.particles[i];
                    glm::vec2 r = (p2.position - p.position);

                    if (dot(r, r) < mPs.supportRadius2) {
                        float r2 = dot(r, r);
                        float ra = sqrt(r2);
                        float h_r = mPs.supportRadius - ra;
                        rhosum += pow((mPs.supportRadius2 - r2), 3.0);
                    }
                }
            }
            p.density = mPs.particleVolume * Lagrangian2dPara::density * mploy * rhosum;
            //std::cout << p.density << std::endl;
            p.pressure = pow((p.density /Lagrangian2dPara::density)-1, Lagrangian2dPara::exponent) * Lagrangian2dPara::density;
        }

        glm::vec2 mul(glm::vec2 a, double b) {
            return glm::vec2(a.x * b, a.y * b);
        }

        void Solver::calx(ParticleInfo2d& p) {
            int bid;
            glm::uvec2 ran;
            glm::vec2 psum(0), vsum(0);
            int sz = int(sqrt(mPs.blockExtens.size()));
            int count = 0;
            float mp = 45.0 / (3.14159f * pow(mPs.supportRadius, 6));
            for (int off : mPs.blockIdOffs) {
                bid = p.blockId + off;
                if (p.blockId % sz == 0 && (bid+1)%sz==0 || (p.blockId + 1) % sz == 0 && bid % sz == 0) {
                    continue;
                }
                ran = mPs.blockExtens[bid];

                if (ran.x == ran.y)continue;
                for (int i = ran.x; i < ran.y; i++) {
                    ParticleInfo2d& p2 = mPs.particles[i];
                    glm::vec2 r = (p2.position - p.position);

                    if (dot(r, r) < mPs.supportRadius2 && &p != &p2) {
                        float ra = sqrt(dot(r,r));
                        float h_r = mPs.supportRadius - ra;
                        
                        psum += mul((p.position - p2.position) ,(p.pressure + p2.pressure) / (2 * p.density * p2.density) * pow(h_r, 2.0) / ra);
                        /*std::cout << mul((p.position - p2.position), (p.pressure + p2.pressure) / (2 * p.density * p2.density) * pow(h_r, 2.0) / ra).x << " " << ra
                            << " " << pow(h_r, 2.0) << " " << (p.pressure + p2.pressure) << " " << (p.position - p2.position).x << " " << psum.x<<std::endl;*/
                        vsum += mul((p2.velocity - p.velocity), h_r/(p.density * p2.density));
                        count++;
                    }
                }
            }
            p.accleration =  mul( (psum + vsum * Lagrangian2dPara::viscosity), mPs.particleVolume * Lagrangian2dPara::density/1.0);
            //std::cout << "2"<<p.accleration.x << std::endl;
            //std::cout <<" " << p.blockId << " "<<count << " " << p.density << " " << p.pressure << " " << p.accleration.x << " " << p.accleration.y << std::endl;

        }
    }
}


