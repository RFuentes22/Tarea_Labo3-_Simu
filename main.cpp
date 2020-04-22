#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "display_tools.h"
#include "tools.h"
#include "sel.h"

int main() {
    vector<Matrix> localKs;
    vector<Vector> localbs;

    //global
    Matrix K;
    Vector b;
    Vector T;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- 1 DIMENSI"<<char(224)<<"N\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m);

    crearSistemasLocales(m,localKs,localbs);
    zeroes(K,m.getSize(NODES));
    zeroes(b,m.getSize(NODES));
    assembly(m,localKs,localbs,K,b);

    applyNeumann(m,b);
    applyDirichlet(m,K,b);
    zeroes(T,b.size());

    calculate(K,b,T);

    cout << "La respuesta es: " << endl;
    showVector(T);


    return 0;
}
