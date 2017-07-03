#include "ecompass.h"

inline void setRotationAngles(uint16_t Gx, uint16_t Gy, uint16_t Gz, rotAngle *protAngle)
{
  protAngle->phi = atan2(Gx, sqrt(Gy*Gy + Gz*Gz));
  protAngle->psi = atan2(Gy, sqrt(Gx*Gx + Gz*Gz));
}

inline void setMagneticFileldVector(uint16_t Bx, uint16_t By, uint16_t Bz, uvector_t *pBp)
{
  pBp->x = Bx;
  pBp->y = By;
  pBp->z = Bz;
}

inline void processTrigonometricFunc(trgRotAngle *ptrgRotAngle, rotAngle *protAngle)
{
  ptrgRotAngle->sin_phi = sin(protAngle->phi);
  ptrgRotAngle->cos_phi = cos(protAngle->phi);
  ptrgRotAngle->sin_theta = sin(protAngle->theta);
  ptrgRotAngle->cos_theta = cos(protAngle->theta);  
}

inline void processAzimuth(uvector_t *pBp, trgRotAngle *ptrgRotAngle, rotAngle *protAngle) 
{
  double Bfx, Bfy;
  
  // equation in matrix form
  //                                    cos(delta)  
  // Bp = Rx(phi)*Ry(theta)*Rz(psi)*B*(     0     ) 
  //                                    sin(delta)
  
  //calculate derotated Bfx
  Bfx = (pBp->x * ptrgRotAngle->cos_theta + pBp->y * ptrgRotAngle->sin_theta
         * ptrgRotAngle->sin_phi + pBp->z * ptrgRotAngle->sin_theta
           * ptrgRotAngle->cos_phi);
  
  //calculate derotated Bfy
  Bfy = (pBp->z * ptrgRotAngle->sin_phi - pBp->y * ptrgRotAngle->cos_phi);
  
  protAngle->theta = atan2(-Bfy, Bfx);
}

