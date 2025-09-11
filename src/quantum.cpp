#include "../includes/quantum.h"
#include <iostream>
#include <cstdlib>
#include <memory>
#include <cmath>

Qdit::Qdit(unsigned dimension) : amplitudes(dimension, 1), compositeSystem(NULL){
    for(int i = 0; i < dimension; i++)
        amplitudes.setValue(
            std::complex<double>(
                std::rand(),std::rand()), i);
    amplitudes.normalize();
}

Qdit::Qdit(ComplexMatrix amplitudes) : amplitudes(amplitudes), compositeSystem(NULL){
    setAmplitudes(amplitudes); 
}

Qdit Qdit::apply(UnitaryOperator U){
    setAmplitudes(U * ket());
    if(compositeSystem)
        compositeSystem->apply(U, compositeSystemIndex);
    return *this;
}

void Qdit::setAmplitudes(ComplexMatrix amplitudes){
    double sum = 0;
    for(int i = 0; i < amplitudes.count(); i++)
        sum += std::pow(std::abs(amplitudes.getValue(i)), 2);
    if(sum - 1.0 > AMPLITUDE_EPSILON){
        std::cout << "ERROR: Invalid qdit amplitudes!" << std::endl;
        exit(0);
    }
    this->amplitudes = amplitudes;
}

void Qdit::setCompositeSystem(CompositeSystem &compositeSystem, int index){
    this->compositeSystem = &compositeSystem;
    compositeSystemIndex = index;
}

ComplexMatrix Qdit::ket(){
    return amplitudes;
}

ComplexMatrix Qdit::bra(){
    return amplitudes.getConjugateTranspose();
}

std::complex<double> Qdit::getAmplitude(unsigned i){
    return amplitudes.getValue(i);
}

ComplexMatrix Qdit::outerProduct(Qdit q1, Qdit q2){
    return q1.ket() * q2.bra(); 
}

std::complex<double> Qdit::braket(Qdit q1, Qdit q2){
    return (q1.bra() * q2.ket()).getValue(0);
}

Qbit::Qbit() : Qdit(2){}
Qbit::Qbit(ComplexMatrix amplitudes) : Qdit(amplitudes){
    if(amplitudes.count() != 2){
        std::cout << "ERROR: Invalid qbit amplitudes!" << std::endl;
        exit(0);
    }    
}
Qbit::Qbit(std::complex<double> alpha, std::complex<double> beta) : Qdit(2){
    if(std::pow(std::abs(alpha), 2) + std::pow(std::abs(beta), 2) - 1.0 > AMPLITUDE_EPSILON){
        std::cout << "ERROR: Invalid qbit amplitudes!" << std::endl;
        exit(0);
    }
    amplitudes.setValue(alpha, 0);
    amplitudes.setValue(beta, 1);
}

void Qbit::colapse(ComplexMatrix amplitudes){
    setAmplitudes(amplitudes);
    if(compositeSystem)
        compositeSystem->colapse();
}

Qbit Qbit::zero()   { return Qbit(1, 0); }
Qbit Qbit::one()    { return Qbit(0, 1); }
Qbit Qbit::plus()   { return Qbit(1/std::sqrt(2), 1/std::sqrt(2)); }
Qbit Qbit::minus()  { return Qbit(1/std::sqrt(2),-1/std::sqrt(2)); }

std::ostream& operator<<(std::ostream& os, const Qbit& q){
    Qbit temp = q;
    std::complex<double> a = temp.getAmplitude(0);
    std::complex<double> b = temp.getAmplitude(1);   
    double x = 2*(a*b).real();
    double y = 2*(a*b).imag();
    double z = std::pow(std::abs(a), 2) - std::pow(std::abs(b), 2);

    char c = 'X';
    char circle[] = {
        ' ', ' ', ' ', '.', ':', '\'', '\'', '\'', ':', '.', ' ', ' ', ' ',
        ' ', '.', '\'', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '`', '.', ' ',
        ':', '\'', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '`', ':',
        ':', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ':',
        ':', '.', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '.', ':',
        ' ', '`', '.', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '.', '\'', ' ',
        ' ', ' ', ' ', '`', ':', '.', '.', '.', ':', '\'', ' ', ' ', ' ', 
    };

    os << "x         |0 >        " << " | ";
    os << "y         |+ >        " << " | ";
    os << "z         |i >        " << std::endl;

    for(int i = 0; i < 7; i++){
        if(i == 3)
            os << "|-i >";
        else
            os << "     ";
        for(int j = 0; j < 13; j++){
            if((3 - i) == std::round(z*3) && (j - 6) == std::round(y*6))
                os << c;
            else
                os << circle[i*13 + j];
        }
        if(i == 3)
            os << "|i >";
        else
            os << "    ";

        os << " | ";

        if(i == 3)
            os << " |1 >";
        else
            os << "     ";
        for(int j = 0; j < 13; j++){
            if((3 - i) == std::round(x*3) && (j - 6) == std::round(z*6))
                os << c;
            else
                os << circle[i*13 + j];
        }
        if(i == 3)
            os << "|0 >";
        else
            os << "    ";

        os << " | ";

        if(i == 3)
            os << " |- >";
        else
            os << "     ";
        for(int j = 0; j < 13; j++){
            if((3 - i) == std::round(y*3) && (j - 6) == std::round(x*6))
                os << c;
            else
                os << circle[i*13 + j];
        }
        if(i == 3)
            os << "|+ >";
        else
            os << "    ";

        os << std::endl;
    }

    os << "          |1 >        " << " | ";
    os << "          |- >        " << " | ";
    os << "         |-i >        " << std::endl;
    os << std::endl;

    return os;
}

UnitaryOperator::UnitaryOperator(ComplexMatrix matrix) : ComplexMatrix(matrix){}

Qdit UnitaryOperator::apply(Qdit &q){
    q.apply((*this));
    return q;
}

double UnitaryOperator::measure(Qbit &q){
    double r = (double)std::rand() / RAND_MAX;
    double p = std::pow(std::abs((getEigenVector(0).getConjugateTranspose()*q.ket()).getValue(0)),2);
    if(p >= r){
        q.colapse(getEigenVector(0));
        return getEigenValue(0);
    }
    q.colapse(getEigenVector(1));
    return getEigenValue(1);
}

double UnitaryOperator::getEigenValue(unsigned i){
    if(std::abs((int)i) > 2){
        std::cout << "ERROR: Eigen value index out off operator bounds!" << std::endl;
        exit(0);
    }
    std::complex<double> p = getValue(0)*getValue(3) - getValue(1)*getValue(2);
    std::complex<double> m = getTrace()/2.0;
    return (m + (1.0 - i*2) * std::sqrt(m*m - p)).real(); 
}

ComplexMatrix UnitaryOperator::getEigenVector(unsigned i){
    ComplexMatrix eigenVector(2,1);
    double lambda = getEigenValue(i);
    std::complex<double> a = getValue(0);
    std::complex<double> b = getValue(1);
    std::complex<double> c = getValue(2);
    std::complex<double> d = getValue(3);     
    if(c != 0.0){
        eigenVector.setValue(lambda - d, 0);
        eigenVector.setValue(c, 1);
    }
    else if(b != 0.0){
        eigenVector.setValue(b, 0);
        eigenVector.setValue(lambda - a, 1);
    }
    else{
        eigenVector.setValue(1.0 - i, 0);
        eigenVector.setValue((double)i, 1);
    }
    return eigenVector.normalize();
}

double UnitaryOperator::expectation(Qbit q){
    return (q.bra() * (*this) * q.ket()).getValue(0).real();
}

UnitaryOperator UnitaryOperator::I(){
    ComplexMatrix mat = Qbit::zero().ket()*Qbit::zero().bra() + 
                        Qbit::one().ket()*Qbit::one().bra();
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::X(){
    ComplexMatrix mat = Qbit::zero().ket()*Qbit::one().bra() + 
                        Qbit::one().ket()*Qbit::zero().bra();
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::Z(){
    ComplexMatrix mat = Qbit::plus().ket()*Qbit::minus().bra() + 
                        Qbit::minus().ket()*Qbit::plus().bra();
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::Y(){
    ComplexMatrix mat = std::complex<double>(0,1) * X() * Z();
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::H(){
    ComplexMatrix mat = Qbit::plus().ket()*Qbit::zero().bra() + 
                        Qbit::minus().ket()*Qbit::one().bra();
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::CNOT(){
    ComplexMatrix mat = ComplexMatrix::tensorProduct(Qbit::zero().ket()*Qbit::zero().bra(), I()) + 
                        ComplexMatrix::tensorProduct(Qbit::one().ket()*Qbit::one().bra(), X());
    return UnitaryOperator(mat);
}

UnitaryOperator UnitaryOperator::Ry(double theta){
    ComplexMatrix mat(2,2);
    mat.setValue(std::cos(theta/2), 0);
    mat.setValue(-std::sin(theta/2),1);
    mat.setValue(std::sin(theta/2), 2);
    mat.setValue(std::cos(theta/2), 3);
    return UnitaryOperator(mat);
}

CompositeSystem::CompositeSystem(Qbit &q1, Qbit &q2) : Qdit(4), q1(&q1), q2(&q2){
    q1.setCompositeSystem(*this, 0);
    q2.setCompositeSystem(*this, 1);
    setAmplitudes(ComplexMatrix::tensorProduct(q1.ket(),q2.ket()));
}

void CompositeSystem::colapse(){
    std::complex<double> a1, a2;
    ComplexMatrix amp1(2,1);
    ComplexMatrix amp2(2,1);
    ComplexMatrix temp1(4,1);
    ComplexMatrix temp2(4,4);
    temp1 = ComplexMatrix::tensorProduct(q1->ket(),q2->ket());

    for(int i = 0; i < 4; i++)
        temp2.setValue(temp1.getValue(i),i,i);
    setAmplitudes((temp2*ket()).normalize());

    if(getAmplitude(0) != 0.0 && getAmplitude(2) != 0.0){
        a1 = getAmplitude(0)/getAmplitude(2);
        amp1.setValue(a1,0);
        amp1.setValue(1,1);
    }
    else if(getAmplitude(1) != 0.0 && getAmplitude(3) != 0.0){
        a1 = getAmplitude(1)/getAmplitude(3);
        amp1.setValue(a1,0);
        amp1.setValue(1,1);
    }
    else if(getAmplitude(0) == 0.0 && getAmplitude(1) == 0.0){
        amp1.setValue(0,0);
        amp1.setValue(1,1);
    }
    else{
        amp1.setValue(1,0);
        amp1.setValue(0,1);
    }
    amp1.normalize();
    q1->setAmplitudes(amp1);

    if(getAmplitude(0) != 0.0 && getAmplitude(1) != 0.0){
        a2 = getAmplitude(0)/getAmplitude(1);
        amp2.setValue(a2,0);
        amp2.setValue(1,1);
    }
    else if(getAmplitude(2) != 0.0 && getAmplitude(3) != 0.0){
        a2 = getAmplitude(2)/getAmplitude(3);
        amp2.setValue(a2,0);
        amp2.setValue(1,1);
    }
    else if(getAmplitude(0) == 0.0 && getAmplitude(2) == 0.0){
        amp2.setValue(0,0);
        amp2.setValue(1,1);
    }
    else{
        amp2.setValue(1,0);
        amp2.setValue(0,1);
    }
    amp2.normalize();
    q2->setAmplitudes(amp2);
}

void CompositeSystem::applyCNOT(){
    ComplexMatrix mat = ComplexMatrix(2,2);
    mat.setValue((q1->bra() * Qbit::zero().ket()).getValue(0),0);
    mat.setValue((q1->bra() * Qbit::one().ket()).getValue(0),1);
    mat.setValue((q1->bra() * Qbit::one().ket()).getValue(0),2);
    mat.setValue((q1->bra() * Qbit::zero().ket()).getValue(0),3);
    q2->setAmplitudes(mat * q2->ket());
    setAmplitudes(UnitaryOperator::CNOT()*ket());
}

void CompositeSystem::apply(UnitaryOperator U, int index){
    if(!index)
        setAmplitudes(ComplexMatrix::tensorProduct(U,UnitaryOperator::I())*ket());
    else
        setAmplitudes(ComplexMatrix::tensorProduct(UnitaryOperator::I(),U)*ket());
}

void CompositeSystem::apply(UnitaryOperator U1, UnitaryOperator U2){
    U1.apply(*q1);
    U2.apply(*q2);
}