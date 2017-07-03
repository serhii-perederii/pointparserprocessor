#ifndef ECOMPASS_H
#define ECOMPASS_H

#include <stdint.h>
#include <math.h>
#include "magnetometer.h"

typedef struct
{
  float theta;
  float phi;
  float psi;
} rotAngle; //rotation angles struct

typedef struct
{
  float cos_theta;
  float sin_theta;
  float cos_phi;
  float sin_phi;
} trgRotAngle; //rotation angles trigonometric function struct 

void setMagneticFileldVector(uint16_t Bx, uint16_t By, uint16_t Bz, uvector_t *pBp);
void setRotationAngles(uint16_t Gx, uint16_t Gy, uint16_t Gz, rotAngle *protAngle);
void processTrigonometricFunc(trgRotAngle *ptrgRotAngle, rotAngle *protAngle);
void processAzimuth(uvector_t *pBp, trgRotAngle *ptrgRotAngle, rotAngle *protAngle);

#endif // ECOMPASS_H