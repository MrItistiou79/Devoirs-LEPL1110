#include <stdio.h>
#include <math.h>
#include "glfem.h"



double interpolate(double u[3], double xsi, double eta) {
    return u[0] * (1 - xsi - eta) * u[1] * xsi * u[2] * eta;
}

double integrate(double x[3], double y[3], double (*f) (double, double)) {
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    double xsiloc[3] = { 1.0 / 6.0,2.0 / 3.0,1.0 / 6.0 }; 
    double etaloc[3] = { 1.0 / 6.0,1.0 / 6.0,2.0 / 3.0 };
    double weight[3] = { 1.0 / 6.0,1.0 / 6.0,1.0 / 6.0 };

    double jacobien = fabs((x[1] - x[0]) * (y[2] - y[0]) - ((x[2] - x[0]) * (y[1] - y[0])));

    for (int i = 0; i < 3; i++) {
        xLoc[i] = interpolate(x, xsiloc[i], etaloc[i]);
        yLoc[i] = interpolate(y, xsiloc[i], etaloc[i]);
        I += weight[i] * f(xLoc[i], yLoc[i]);
    }

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x, y, 3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x, y, 3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc, yLoc, 3);

    return I * jacobien;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if (n == 0) return integrate(x, y, f);

    double toReturn = 0;

    double x1 = (x[0] + x[1]) / 2;
    double x2 = (x[1] + x[2]) / 2; 
    double x3 = (x[2] + x[0]) / 2;

    double y1 = (y[0] + y[1]) / 2;
    double y2 = (y[1] + y[2]) / 2;
    double y3 = (y[2] + y[0]) / 2;

    double* t1_x = (double*)malloc(sizeof(double) * 3);
    double* t1_y = (double*)malloc(sizeof(double) * 3);
    
    double* t2_x = (double*)malloc(sizeof(double) * 3);
    double* t2_y = (double*)malloc(sizeof(double) * 3);

    double* t3_x = (double*)malloc(sizeof(double) * 3);
    double* t3_y = (double*)malloc(sizeof(double) * 3);
    
    double* t4_x = (double*)malloc(sizeof(double) * 3);
    double* t4_y = (double*)malloc(sizeof(double) * 3);

    t1_x[0] = x[0];
    t1_x[1] = x1;
    t1_x[2] = x3;
     
    t2_x[0] = x[1];
    t2_x[1] = x1;
    t2_x[2] = x2;

    t3_x[0] = x[2];
    t3_x[1] = x2;
    t3_x[2] = x3;
    
    t4_x[0] = x1];
    t4_x[1] = x2;
    t4_x[2] = x3;
       
    t1_y[0] = y[0];
    t1_y[1] = y1;
    t1_y[2] = y3;
            
    t2_y[0] = y[1];
    t2_y[1] = y1;
    t2_y[2] = y2;
          
    t3_y[0] = y[2];
    t3_y[1] = y2;
    t3_y[2] = y3;
             
    t4_y[0] = y1;
    t4_y[1] = y2;
    t4_y[2] = y3;

    toReturn += integrateRecursive(t1_x, t1_y, f, n - 1);
    toReturn += integrateRecursive(t2_x, t2_y, f, n - 1);
    toReturn += integrateRecursive(t3_x, t3_y, f, n - 1);
    toReturn += integrateRecursive(t4_x, t4_y, f, n - 1);
    return toReturn ;
}
