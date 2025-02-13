#include "fluid3d/Lagrangian/include/Solver.h"

namespace FluidSimulation
{

	namespace Lagrangian3d
	{
		Solver::Solver(ParticleSystem3d &ps) : mPs(ps)
		{
		}

		void Solver::solve()
		{
			// TODO
			// Solves the fluid simulation by performing some steps, which may include:
			// 1. compute density 
			// 2. compute press
			// 3. compute accleration
			// 4. update velocity and position
			// 5. check boundary
			// 6. update block id
			// ...

			calpre();
			calacc();
			calx();
			mPs.updateBlockInfo();
		}

		void Solver::calpre() {
			double poly = 315.0f / (64.0f * 3.141592f * pow(mPs.supportRadius, 9));
			for (auto& particle : mPs.particles) {
				float sum1 = 0.0;

				int blockId = particle.blockId;
				for (int offset : mPs.blockIdOffs) {
					int neighborBlockId = blockId + offset;
					if (neighborBlockId >= 0 && neighborBlockId < mPs.blockExtens.size()) {
						glm::uvec2 tmp = mPs.blockExtens[neighborBlockId]; // 假设blockExtens现在是三维的
						for (int i = tmp[0]; i < tmp[1]; ++i) {
							particle3d& neighbor = mPs.particles[i];
							if (&neighbor != &particle) {
								glm::vec3 r = particle.position - neighbor.position;
								float dist2 = glm::dot(r, r);
								if (dist2 < mPs.supportRadius2) {
									sum1 += pow((mPs.supportRadius2 - dist2), 3.0);
									//std::cout << sum1 << std::endl;
								}
							}
						}
					}
				}

				particle.density = poly * Lagrangian3dPara::density * mPs.particleVolume * 3.141592f * sum1;
				particle.pressure = 1.0f * (particle.density - Lagrangian3dPara::density);
			}
		}

		void Solver::calacc() {
			double spicky = -45.0f / (3.141592f * pow(mPs.supportRadius, 6));
			double vis = 45.0f / (3.141592f * pow(mPs.supportRadius, 6));
			for (auto& particle : mPs.particles) {
				glm::vec3 accel_sum = { 0.0, 0.0, 0.0 };

				int blockId = particle.blockId;
				for (int offset : mPs.blockIdOffs) {
					int neighborBlockId = blockId + offset;
					if (neighborBlockId >= 0 && neighborBlockId < mPs.blockExtens.size()) {
						glm::uvec2 tmp = mPs.blockExtens[neighborBlockId]; // 假设blockExtens现在是三维的
						for (int i = tmp[0]; i < tmp[1]; ++i) {
							particle3d& neighbor = mPs.particles[i];
							if (&neighbor != &particle) {
								glm::vec3 r = particle.position - neighbor.position;
								float rr = glm::length(r);
								if (rr > 1e-6 && rr < mPs.supportRadius) { // 避免除零并限制在支持半径内
									float h_r = mPs.supportRadius - rr;
									glm::vec3 vi_vj = particle.velocity - neighbor.velocity;

									float pterm = spicky * h_r * h_r * (particle.pressure + neighbor.pressure) / (2.0f * particle.density * neighbor.density);
									accel_sum -= r * pterm / rr;

									float vterm = vis * Lagrangian3dPara::viscosity * h_r / (particle.density * neighbor.density);
									accel_sum += vi_vj * vterm;
								}
							}
						}
					}
				}
				particle.accleration = -glm::vec3(Lagrangian3dPara::gravityX, Lagrangian3dPara::gravityY, Lagrangian3dPara::gravityZ) + accel_sum * Lagrangian3dPara::density * mPs.particleVolume * 3.141592f;
			}
		}

		void Solver::calx() {
			int i = 0;
			//std::cout << mPs.m_kernelSpiky <<" "<< mPs.m_kernelPoly6 <<" "<< mPs.m_kernelViscosity << std::endl;
			for (auto& particle : mPs.particles) {
				glm::vec3 oldPos = particle.position;
				//std::cout << "(" << particle.accleration[0] << "," << particle.accleration[1] << "," << particle.accleration[2] << ") ";

				particle.velocity += Lagrangian3dPara::dt * particle.accleration;
				particle.position += Lagrangian3dPara::dt * particle.velocity;
				//std::cout << "[" << particle.position[0] << "," << particle.position[1] << "," << particle.position[2] << "] "<<std::endl;


				if (glm::dot(oldPos - particle.position, oldPos - particle.position) > 0.00001) {
					//std::cout << i << ": " << particle.acceleration.x << " " << particle.acceleration.y << " " << particle.acceleration.z << " " 
					//          << particle.density << " " << particle.pressure << std::endl;
				}
				i++;

				// 边界条件处理
				if (particle.position.x < mPs.lowerBound.x) {
					particle.position.x = mPs.lowerBound.x;
					particle.velocity.x *= -Lagrangian3dPara::velocityAttenuation;
				}
				if (particle.position.y < mPs.lowerBound.y) {
					particle.position.y = mPs.lowerBound.y;
					particle.velocity.y *= -Lagrangian3dPara::velocityAttenuation;
				}
				if (particle.position.z < mPs.lowerBound.z) {
					particle.position.z = mPs.lowerBound.z;
					particle.velocity.z *= -Lagrangian3dPara::velocityAttenuation;
				}
				if (particle.position.x > mPs.upperBound.x) {
					particle.position.x = mPs.upperBound.x;
					particle.velocity.x *= -Lagrangian3dPara::velocityAttenuation;
				}
				if (particle.position.y > mPs.upperBound.y) {
					particle.position.y = mPs.upperBound.y;
					particle.velocity.y *= -Lagrangian3dPara::velocityAttenuation;
				}
				if (particle.position.z > mPs.upperBound.z) {
					particle.position.z = mPs.upperBound.z;
					particle.velocity.z *= -Lagrangian3dPara::velocityAttenuation;
				}
			}
		}
	}
}