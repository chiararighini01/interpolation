#include "Interpolation/generic_grid.hh"

namespace Interpolation
{
namespace Generic
{
    //costruttore a partire da un vettore di punti
StandardGrid::StandardGrid(const vector_d &input) : _tj(input) //after: inizializer list (copy of input vector, defines the _tj before entering in the main of the constructor)
{
    if(input.size()<=1){
        throw std::invalid_argument("[Generic::StandardGrid]: input vector cannot be less than two elements.");
    }
    
    std::sort(_tj.begin(), _tj.end()); //riordina gli elementi di _tj (lo cambia proprio)

    //begin e end sono iteratori, front e back sono elementi

    //the input vector ha to contain points in [-1,1]
    if (_tj.front() < -1.) { 
      throw std::runtime_error("Invalid vector");
    }

    if (_tj.back() > +1.) {
        throw std::runtime_error("Invalid vector");
    }

    //the first point has to be -1 and the last one 1
    if (std::abs(_tj.front() - (-1.)) > 1.0e-12) {  //se il primo punto non è -1
        _tj.insert(_tj.begin(), -1.);               //aggiungi il punto -1 all'inizio
    } else {                                        //se il primo punto è circa -1
        _tj.front() = -1.;                          //fallo diventare esattamente -1
    }
    if (std::abs(_tj.back() - 1.) > 1.0e-12) {      //lo stesso con il +1 finale
        _tj.push_back(+1.);
    } else {
        _tj.back() = +1.;
    }

    _p=_tj.size()-1;
    _lambdaj.resize(_p+1,1.);                        //inizializza i pesi lambda j a un vettore lungo p+1 di 1

    //lambda_j=prod for i!=j of tj-ti
    for (size_t j = 0; j <= _p; j++) {
      for (size_t i = 0; i < j; i++) {              //invece di fare i 2 loop potrei farne uno solo con la condizione che i sia diverso da j
         _lambdaj[j] *= _tj[j] - _tj[i];
      }
      for (size_t i = j + 1; i <= _p; i++) {
         _lambdaj[j] *= _tj[j] - _tj[i];
      }
      _lambdaj[j] = 1. / _lambdaj[j];
   }

    //Derivative matrix
    _Dij.resize(_p + 1, std::vector<double>(_p+1,0.)); //inizializza un vettore di vettori pieno di 0

    for(size_t i=0;i<=_p;i++){
        //Diagonal elements
        for(size_t n = 0; n<i; n++){
            _Dij[i][i] += 1./ (_tj[i]-_tj[n]);
        }
        for(size_t n = i+1; n<=_p; n++){
            _Dij[i][i] += 1./ (_tj[i]-_tj[n]);
        }
        //Off-diagonal elements
        for(size_t j=0; j<=_p;j++){
            if(i==j) continue;
            _Dij[i][j]=1./(_tj[i]-_tj[j]);
            for(size_t k=0; k<=_p;k++){
                if(k==i || k==j) continue;
                _Dij[i][j]*=(_tj[j]-_tj[k])/(_tj[i]-_tj[k]);
            }

        }
    }
}

//costruttore a partire da una funzione a 2 variabili e dal numero di punti ?? (cos'è questa fne a 2 variabili? se ne usa solo una tra l'altro...)
StandardGrid::StandardGrid(const std::function<double(size_t, size_t)> &fnc, size_t p): _p(p)
{
    //Input checks
    if(std::fabs(fnc(0,p) - (-1.)) > 1.0e-12){
        throw std::invalid_argument("[Generic::StandardGrid] lower bound of the grid must be -1");

    }

    //Weights
    _lambdaj.resize(_p+1,1.);                        //inizializza i pesi lambda j a un vettore lungo p+1 di 1

    for (size_t j = 0; j <= _p; j++) {
      for (size_t i = 0; i < j; i++) {              //invece di fare i 2 loop potrei farne uno solo con la condizione che i sia diverso da j
         _lambdaj[j] *= _tj[j] - _tj[i];
      }
      for (size_t i = j + 1; i <= _p; i++) {
         _lambdaj[j] *= _tj[j] - _tj[i];
      }
      _lambdaj[j] = 1. / _lambdaj[j];
    }

    //Derivative matrix
    _Dij.resize(_p + 1, std::vector<double>(_p+1,0.)); //inizializza un vettore di vettori pieno di 0

    for(size_t i=0;i<=_p;i++){
        //Diagonal elements
        for(size_t n = 0; n<i; n++){
            _Dij[i][i] += 1./ (_tj[i]-_tj[n]);
        }
        for(size_t n = i+1; n<=_p; n++){
            _Dij[i][i] += 1./ (_tj[i]-_tj[n]);
        }
        //Off-diagonal elements
        for(size_t j=0; j<=_p;j++){
            if(i==j) continue;
            _Dij[i][j]=1./(_tj[i]-_tj[j]);
            for(size_t k=0; k<=_p;k++){
                if(k==i || k==j) continue;
                _Dij[i][j]*=(_tj[j]-_tj[k])/(_tj[i]-_tj[k]);
            }

        }
    }

}
double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end, STRATEGY str) const
{
    //Input checks
    if(t<-1||t>1 ){
        throw std::invalid_argument("[Generic::StandardGrid::interpolate] t = "+ std::to_string(t)+ " not in [-1,1]");
    }
    if((end-start)!=_p){
        throw std::invalid_argument("[Generic::StandardGrid::interpolate] fj of the wrong size: ["+ std::to_string(start)+ " ,"+std::to_string(end)+"]");
    }

    //strategy selection
    switch (str) {
        case STRATEGY::NAIVE: {
            double sum = 0;
            for (size_t i = start, j = 0; i <= end; i++, j++) {  //i va da start a end mentre j va contemporaneamente da 0 a p
                sum += poli_weight(t, j) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::FBF: {
            double monic = 1.;
            for (size_t i = 0; i <= _p; i++) {
                monic *= (t - _tj[i]);
            }
            double sum = 0.;
            for (size_t i = start, j = 0; i <= end; i++, j++) {
                sum += poli_weight_fbf(t, j, monic) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::SBF: {
            double den = 0;
            for (size_t l = 0; l <= _p; l++) {
                if (fabs(t - _tj[l]) < 1.0e-15) return fj[l + start]; //se valuto esattamente in un punto della griglia
                den += _lambdaj[l] / (t - _tj[l]);
            }

            double sum = 0;
            for (size_t i = start, j = 0; i <= end; i++, j++) {
                sum += poli_weight_sbf(t, j, den) * fj[i];
            }
            return sum;
        }
    }
}

double StandardGrid::poli_weight(double t, size_t j) const
{
    //input checks 
    if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare

    double num=0.;
	for (size_t i=0; i<=_p; i++){
        if(i==j) continue;
		if (std::abs(t-_tj[i]) < 1.0e-15 ) return 0.; //if the difference is really small
		num*=(t-_tj[i]);
	}
	
	double res=0.;
	res = _lambdaj[j]*num;
	return res;
}

double StandardGrid::poli_weight_fbf(double t, size_t j, double monic) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare
    return monic*_lambdaj[j]/(t-_tj[j]);
}

double StandardGrid::poli_weight_fbf(double t, size_t j) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare
    double monic = 1.;
    for (size_t i = 0; i <= _p; i++) {
        monic *= (t - _tj[i]);
    }
       
    return monic*_lambdaj[j]/(t-_tj[j]);
}

double StandardGrid::poli_weight_sbf(double t, size_t j, double den) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare
    return _lambdaj[j]/den/(t-_tj[j]);
}

double StandardGrid::poli_weight_sbf(double t, size_t j) const
{
    if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare
    
    double den = 0;
    for (size_t l = 0; l <= _p; l++) {
        if (fabs(t - _tj[l]) < 1.0e-15) return 0.; //se valuto esattamente in un punto della griglia
        den += _lambdaj[l] / (t - _tj[l]);
    }
    
    return _lambdaj[j]/den/(t-_tj[j]);
}

vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
    vector_d fj(_p+1, 0.);;
	for(size_t i=0; i<=_p; i++){
		fj[i] = fnc(_tj[i]);
	}
	return fj;
}

double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end,
                          STRATEGY str) const
{
    //Input checks
    if(t<-1||t>1 ){
        throw std::invalid_argument("[Generic::StandardGrid::interpolate_der] t = "+ std::to_string(t)+ " not in [-1,1]");
    }
    if((end-start)!=_p){
        throw std::invalid_argument("[Generic::StandardGrid::interpolate_der] fj of the wrong size: ["+ std::to_string(start)+ " ,"+std::to_string(end)+"]");
    }
    //strategy selection
    switch (str) {
        case STRATEGY::NAIVE: {
            double sum = 0;
            for (size_t i = start, j = 0; i <= end; i++, j++) {  //i va da start a end mentre j va contemporaneamente da 0 a p
                sum += poli_weight_der(t, j) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::FBF: {
            double monic = 1.;
            for (size_t i = 0; i <= _p; i++) {
                monic *= (t - _tj[i]);
            }
            double sum = 0.;
            for (size_t i = start, j = 0; i <= end; i++, j++) {
                sum += poli_weight_fbf_der(t, j, monic) * fj[i];
            }
            return sum;
            break;
        }
        case STRATEGY::SBF: {
            double den = 0;
            for (size_t l = 0; l <= _p; l++) {
                if (fabs(t - _tj[l]) < 1.0e-15) return fj[l + start]; //se valuto esattamente in un punto della griglia
                den += _lambdaj[l] / (t - _tj[l]);
            }

            double sum = 0;
            for (size_t i = start, j = 0; i <= end; i++, j++) {
                sum += poli_weight_sbf_der(t, j, den) * fj[i];
            }
            return sum;
        }
    }

}

double StandardGrid::interpolate_der_v2(double t, const vector_d &fj, size_t start, size_t end, STRATEGY str) const
{
    vector_d ftilde = fj; //creo una copia per modificarla
    apply_D(ftilde, start, end);
    return interpolate(t, ftilde, start, end, str);
}


double StandardGrid::poli_weight_der(double t, size_t j) const
{
	//input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight_der]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }

	double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight(t, i);
	}
	return sum;
}

double StandardGrid::poli_weight_fbf_der(double t, size_t j) const
{
    //input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight_fbf_der]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
    }

    double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight_fbf(t, i);
	}
	return sum;
}

double StandardGrid::poli_weight_fbf_der(double t, size_t j, double monic) const
{
    //input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight_fbf_der]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
    }

    double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight_fbf(t, i, monic);
	}
	return sum;
}

double StandardGrid::poli_weight_sbf_der(double t, size_t j) const
{
    //input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
    }

    double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight_sbf(t, i);
	}
	return sum;
}

/// Same, for second barycentric formula with precomputed denominator
double StandardGrid::poli_weight_sbf_der(double t, size_t j, double den) const
{
    //input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight_sbf_der]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
    }

    double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight_sbf(t, i, den);
	}
	return sum;
}

void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const
{
	//input checks
   if (end - start != _p) {
      throw std::invalid_argument("[StandardGrid::apply_D]: cannot apply "
                                  "derivative matrix to partial vector.");
   }

   vector_d temp((end - start + 1), 0.);
   for (size_t i = 0; i <= _p; i++) {
      for (size_t j = 0, k = start; k <= end; k++, j++) {
         temp[i] += fj[k] * _Dij[j][i];
      }
   }

   for (size_t i = start; i <= end; i++) {
      fj[i] = temp[i - start];
   }
}


} // namespace Generic
} // namespace Interpolation