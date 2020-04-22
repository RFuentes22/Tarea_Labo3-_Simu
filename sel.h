
Matrix createLocalK(int element,mesh &m){
    Matrix K;
    Vector row1,row2;

    float e = m.getParameter(CTE_E);
    float a = m.getParameter(CTE_A);
    float l = m.getNode(element + 1).getX() - m.getNode(element).getX();

    row1.push_back((e*a/l)/2); row1.push_back(-(e*a/l)/2);
    row2.push_back(-(e*a/l)/2); row2.push_back((e*a/l)/2);

    K.push_back(row1);
    K.push_back(row2);
   // showMatrix(K);
    return K;
}

Vector createLocalb(int element,mesh &m){
    Vector b;
    float f = m.getParameter(CTE_F);
    float l = m.getNode(element + 1).getX() - m.getNode(element).getX();

    b.push_back((f*l)/2);
    b.push_back((f*l)/2);
    //showVector(b);
    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localks,vector<Vector> &localbs){

    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localks.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localk,Matrix &K){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    K.at(index1).at(index1) += localk.at(0).at(0);
    K.at(index1).at(index2) += localk.at(0).at(1);
    K.at(index2).at(index1) += localk.at(1).at(0);
    K.at(index2).at(index2) += localk.at(1).at(1);
    //showMatrix(K);
}

void assemblyb(element e,Vector localb,Vector &b,mesh &m){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    b.at(index1) += -1*localb.at(0);
    b.at(index2) += -1*localb.at(1);
    //showVector(b);
}

void assembly(mesh &m,vector<Matrix> &localks,Matrix &localbs,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localks.at(i),K);
        assemblyb(e,localbs.at(i),b,m);
    }
}

void applyNeumann(mesh &m,Vector &b){
    float e = m.getParameter(CTE_E);
    float a = m.getParameter(CTE_A);
    for(int i=0;i<m.getSize(NEUMANN);i++){
        condition c = m.getCondition(i,NEUMANN);
        b.at(c.getNode1()-1) += c.getValue();
    }
    //showVector(b);

}

void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    float a = m.getParameter(CTE_A);
    float e = m.getParameter(CTE_E);

    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);
        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}

void calculate(Matrix &K,Vector &b,Vector &T){
    Matrix Kinv;
    inverseMatrix(K,Kinv);
    productMatrixVector(Kinv,b,T);
}