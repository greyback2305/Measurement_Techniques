#include <Wire.h>
#include <SPI.h>
#include "SparkFunLIS3DH.h"
#include <math.h>

LIS3DH myIMU(I2C_MODE, 0x18);

float ax, ay, az;
float aAbs;
float roll, pitch;

void setup() {
  Serial.begin(115200);
  delay(1000);

  Serial.println("Demarrage...");

  if (myIMU.begin() != 0) {
    Serial.println("Erreur : LIS3DH non detecte.");
    while (1);
  }

  Serial.println("LIS3DH detecte.");
  Serial.println("Temps Ax Ay Az Aabs Roll Pitch");
}

void loop() {
  float temps = millis() / 1000.0;

  ax = myIMU.readFloatAccelX() * 1000;
  ay = myIMU.readFloatAccelY() * 1000;
  az = (myIMU.readFloatAccelZ()) * 1000;

  aAbs = sqrt(ax * ax + ay * ay + az * az);

  // Tilt angles in degrees (uses raw az without the +1g offset)
  float az_raw = myIMU.readFloatAccelZ();
  float ax_raw = myIMU.readFloatAccelX();
  float ay_raw = myIMU.readFloatAccelY();

  roll  = atan2(ay_raw, az_raw) * 180.0 / M_PI;
  pitch = atan2(-ax_raw, sqrt(ay_raw * ay_raw + az_raw * az_raw)) * 180.0 / M_PI;

    // In loop(), replace Serial prints with:
  Serial.print(temps); Serial.print(",");
  Serial.print(ax, 3); Serial.print(",");
  Serial.print(ay, 3); Serial.print(",");
  Serial.print(az, 3); Serial.print(",");
  Serial.print(aAbs, 3); Serial.print(",");
  Serial.print(roll, 2); Serial.print(",");
  Serial.println(pitch, 2);

  delay(20);
}