#include <Arduino.h>
#include "AS5600.h"
#include <Wire.h>

AS5600 as5600;

// Timing constants
const int readInterval = 10;                           // Time between readings (ms)
const int averageInterval = 100;                       // Time to calculate average (ms)
const int numSamples = averageInterval / readInterval; // Number of samples in each averaging period

// Variables for calibration and averaging
unsigned long lastReadTime = 0;
unsigned long lastAverageTime = 0;
int angleSum = 0;
int sampleCount = 0;
int calAngle1 = 0, calAngle2 = 0; // Calibration points
float zeroAngle = 0;              // New zero-point in degrees

// Function to calibrate the encoder based on two known angles
void calibrate()
{
  Serial.println("Calibration: Rotate to first known angle and type 'y' to capture.");
  while (Serial.read() != 'y')
  {
  } // Wait for 'y' input
  calAngle1 = as5600.rawAngle();
  Serial.print("Captured first angle: ");
  Serial.println(calAngle1);

  delay(500); // Brief pause for user to adjust position

  Serial.println("Calibration: Rotate to second known angle and type 'y' to capture.");
  while (Serial.read() != 'y')
  {
  } // Wait for 'y' input
  calAngle2 = as5600.rawAngle();
  Serial.print("Captured second angle: ");
  Serial.println(calAngle2);

  // Calculate the midpoint between the two calibration points
  zeroAngle = ((calAngle1 + calAngle2) / 2.0) * (360.0 / 4096.0); // Convert to degrees
  Serial.print("Calibration complete. New zero-point angle: ");
  Serial.println(zeroAngle);
}

void setup()
{
  Serial.begin(9600);
  while (!Serial)
  {
  };

  Wire.begin();

  Serial.print("Beginning: ");
  Serial.println(as5600.begin(1)); // Set direction pin.

  int connected = as5600.isConnected();
  Serial.print("Connected: ");
  Serial.println(connected);

  // Prompt user to calibrate the encoder
  calibrate();
  delay(1000);
}

void loop()
{
  unsigned long currentTime = millis();

  // Read every 10 milliseconds
  if (currentTime - lastReadTime >= readInterval)
  {
    lastReadTime = currentTime;

    // Read raw angle and convert to degrees
    int rawAngle = as5600.rawAngle();
    float angleDegrees = (rawAngle / 4096.0) * 360.0;

    // Adjust angle based on calibration zero-point
    float adjustedAngle = angleDegrees - zeroAngle;
    /*if (adjustedAngle < 0)
    {
      adjustedAngle += 360.0; // Wrap angle to be within 0-360
    }*/

    // Accumulate adjusted angle for averaging
    angleSum += adjustedAngle;
    sampleCount++;
  }

  // Calculate and report average every 100 milliseconds
  if (currentTime - lastAverageTime >= averageInterval)
  {
    lastAverageTime = currentTime;

    // Calculate average angle
    float averageAngle = (sampleCount > 0) ? (float(angleSum) / sampleCount) : 0;

    // Print the average adjusted angle
    Serial.print("Average Adjusted Angle (degrees): ");
    Serial.println(averageAngle);

    // Reset sum and sample count for the next interval
    angleSum = 0;
    sampleCount = 0;
  }
}