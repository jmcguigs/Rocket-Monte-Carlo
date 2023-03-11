// cppStatisticalComputing.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include "..\VectorAlgebra.h"
using namespace std;

using namespace VectorAlgebra;


// BEGIN MISC FUNCTIONS

// get atmospheric density using height above sea level
float getAtmosphericDensity(float HSL)
{
    float referencePressure = 1013250.00;   // reference pressure in Pa
    float referenceTemperature = 288.15;    // reference temperature in K
    float tempLapseRate = -0.0065;          // temperature lapse rate in K/m
    float baseHeight = 0;                   // height of reference level
    float gasConstant = 8.3144598;          // universal gas constant in J/(mol*K)
    float gravAcceleration = 9.80665;       // gravitational acceleration in m/s^2
    float airMolarMass = 0.0289644;         // molar mass of air in kg/mol

    float pressure = referencePressure * pow((referenceTemperature + (HSL - baseHeight) * tempLapseRate) / referenceTemperature, -gravAcceleration * airMolarMass / (gasConstant * tempLapseRate));
    return pressure * (airMolarMass / (gasConstant * referenceTemperature));
}

/*
float polynomialEvaluate(float * coefficients, float x)
{
    float result = 0;
    for (int i = 0; !isnan(coefficients[i]); i++)
    {
        result += coefficients[i] * pow(x, i);
    }

    return result;
}
*/

float polynomialEvaluate(float* coefficients, int order, float x)
{
    float result = 0;
    for (int i = 0; i <= order; i++)
    {
        result += coefficients[i] * pow(x, order - i);
    }

    return result;
}

// END MISC FUNCITONS

class Rocket
{
public:
    Vector3D position = Vector3D(0, 0, 0);
    Vector3D velocity = Vector3D(0, 0, 0);
    Vector3D acceleration = Vector3D(0, 0, -9.81);
    Vector3D orientation = Vector3D(0, 0, 1);

    float radius;                   // rocket tube diameter
    float dragCoefficient;          // drag coefficient of rocket body
    float centerOfPressureOffset;   // distance from center of mass to center of pressure
    float mass;                     // rocket wet mass
    float * thrustCurveCoefficients; // polynomial coefficients of thrust curve
    int thrustCurvePolyOrder;       // order of thrust curve polynomial
    float burnDuration;


    // assign a new position to the rocket
    void setPosition(float x, float y, float z)
    {
        position = Vector3D(x, y, z);
    }

    // set rocket's velocity
    void setVelocity(float x, float y, float z)
    {
        velocity = Vector3D(x, y, z);
    }

    // set rocket's acceleration
    void setAcceleration(float x, float y, float z)
    {
        acceleration = Vector3D(x, y, z);
    }

    // set orientation vector using polar coordinate angles
    void setOrientation(float phi, float theta)
    {
        float p = 3.141592 * phi / 180;
        float t = 3.141592 * theta / 180;

        orientation = Vector3D(sin(t) * cos(p), sin(t) * sin(p), cos(t));
    }


    float getDragForce(float speed)
    {
        //float speed = sqrt(velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z);
        return 0.5 * getAtmosphericDensity(position.z) * speed * speed * dragCoefficient * 3.141592 * radius * radius;
    }

    // simulate rocket launch and find apogee
    float getApogeeAltitude(void)
    {
        float apogee = 0;
        int step = 0;
        float stepSize = 0.001;
        Vector3D thrust = Vector3D(0, 0, 0);
        while (velocity.z >= 0)
        {
            float time = step * stepSize;
            if (time > burnDuration)
            {
                thrust = Vector3D(0, 0, 0);
            }
            else
            {
                //thrust = polynomialEvaluate(thrustCurveCoefficients, time);
                //std::cout << thrust << "\n";
                thrust = orientation * polynomialEvaluate(thrustCurveCoefficients, thrustCurvePolyOrder, time);
                //thrust = orientation * 1500;
            }
            acceleration = Vector3D(0, 0, -9.81) + thrust / mass;
            acceleration = acceleration - Vector3D(0, 0, getDragForce(velocity.z));

            velocity = velocity + acceleration * stepSize;
            //std::cout << velocity.z << "\n";
            position = position + velocity * stepSize;
            apogee = position.z;
            step += 1;
            //std::cout << apogee << "\n";
        }

        return apogee;
    }

    // simulate rocket launch and find apogee
    Vector3D* getFlightPath(void)
    {
        float apogee = 0;
        int step = 0;
        float stepSize = 0.001;
        Vector3D thrust = Vector3D(0, 0, 0);
        while (velocity.z >= 0)
        {
            float time = step * stepSize;
            if (time > burnDuration)
            {
                thrust = Vector3D(0, 0, 0);
            }
            else
            {
                //thrust = polynomialEvaluate(thrustCurveCoefficients, time);
                //std::cout << thrust << "\n";
                thrust = orientation * polynomialEvaluate(thrustCurveCoefficients, thrustCurvePolyOrder, time);
                //thrust = orientation * 1500;
            }
            acceleration = Vector3D(0, 0, -9.81) + thrust / mass;
            acceleration = acceleration - Vector3D(0, 0, getDragForce(velocity.z));

            velocity = velocity + acceleration * stepSize;
            //std::cout << velocity.z << "\n";
            position = position + velocity * stepSize;
            apogee = position.z;
            step += 1;
            //std::cout << apogee << "\n";
        }

        return apogee;
    }

};

int main()
{
    Rocket run1 = Rocket();
    run1.mass = 30;
    run1.setVelocity(0, 0, 5);
    //run1.setOrientation(0, 10);
    run1.radius = 0.01;
    run1.dragCoefficient = 0.1;
    run1.burnDuration = 7;
    float coeffs[5] = { 1.83838, -23.8556, -5.28313, 306.234, 1658.34 };
    //float coeffs[5] = { 1658.34, 306.234, -5.28313, -23.8556, 1.83838 };
    //float coeffs[7] = { 1507.04, 952.044, -640.282, 210.39, -33.2858, 1.57273, 0.0379133 };
    run1.thrustCurveCoefficients = coeffs;
    run1.thrustCurvePolyOrder = 4;

    ofstream outputFile;
    outputFile.open("C:\\Users\\jemcg\\Desktop\\out.txt");

    /*
    // get apogee for a linear sweep of launch angles
    for (int angle = 0; angle <= 60; angle++)
    {
        run1.setOrientation(0, angle);
        outputFile << angle << "," << run1.getApogeeAltitude() << "\n";
        run1.setVelocity(0, 0, 5);
        run1.setAcceleration(0, 0, -9.81);
        run1.setPosition(0, 0, 0);
    }

    */

    // random number generator setup
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    // distribution of launch angles
    std::normal_distribution<> launchAngles{ 0, 2 };
    std::normal_distribution<> wetMass{ 30, 1 };

    for (int n = 0; n < 10000; n++)
    {
        float angle = launchAngles(gen);
        run1.setOrientation(0, angle);
        //run1.mass = wetMass(gen);
        outputFile << angle << "," << run1.getApogeeAltitude() << "\n";
        run1.setVelocity(0, 0, 5);
        run1.setAcceleration(0, 0, -9.81);
        run1.setPosition(0, 0, 0);
    }

    outputFile.close();

    //std::cout << run1.getApogeeAltitude();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
