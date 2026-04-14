#include "Interpolation/chebyshev_grid.hh"
#include <stdexcept>


namespace Interpolation
{
namespace Chebyshev
{
StandardGrid::StandardGrid(size_t p)
{
   _p = p;

   _tj.resize(_p + 1, 0.);
   _betaj.resize(_p + 1, 0.);
   for (size_t i = 0; i <= _p; i++) {
      _tj[i]    = cos(i * M_PI / (static_cast<double>(_p)));
      _betaj[i] = (i % 2 == 0) ? +1 : -1;
      if (i == 0 || i == p) _betaj[i] *= 0.5;
   }

   _Dij.resize(_p + 1, vector_d(_p + 1, 0.));
   _Dij[0][0]   = (2 * _p * _p + 1) / 6.0;
   _Dij[_p][_p] = -_Dij[0][0];
   //Diagonal elements
   for (size_t j = 1; j < _p; j++) {
      _Dij[j][j] = -_tj[j] * 0.5 / (1.0 - _tj[j] * _tj[j]);
   }
   //Not diagonal elements
   for (size_t j = 0; j <= _p; j++) {
      for (size_t k = 0; k <= _p; k++) {
         if (k == j) continue;
         _Dij[j][k] = -(_betaj[k] / _betaj[j]) / (_tj[j] - _tj[k]);
      }
   }
}
   
double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const  
//const non verrà modificato nella funzione, & ha significati diversi in C e in C++, in questo caso significa reference, non copia ma punta all'indirizzo, meglio mettere const davanti per evitare di modificare il vettore originale
{
	//input checks
	if(t<-1 || t>1){
		throw std::domain_error("StandardGrid::interpolate t must be in [-1,1]");
	}
	if(end - start != _p){
		throw std::domain_error("StandardGrid::interpolate end-start should be = to p");
	}

	//precompilation of the denominator
	double den=0;
	for (size_t j=0; j<=_p; j++){
//		if (t==tj[j]) return fj[j+start];
		if (std::abs(t-_tj[j]) < 1.0e-15 ) return fj[j+start]; //if the difference is really small

		den+=_betaj[j]/(t-_tj[j]);
	}
	
	double res= 0.;
	for (size_t i=0; i<=_p; i++){
		res+= poli_weight(t,i, den)* fj[i + start];
	}
	return res;
}

double StandardGrid::poli_weight(double t, size_t j, double den) const
{
	if (std::abs(t-_tj[j]) < 1.0e-15 ) return 1; //se chiedo il peso di un punto della griglia, se lo chiamo da interpolate ho già fatto il check, ma se lo uso da altre parti lo devo ricontrollare

	double res=0.;
	res = _betaj[j]/(t-_tj[j]) / den;
	return res;
}

double StandardGrid::poli_weight(double t, size_t j) const //nelle note sarebbe lj(t)
{
//prima calcolo il denominatore
	double den=0;
	for (size_t j=0; j<=_p; j++){
//		if (t==tj[j]) return fj[j+start];
		if (std::abs(t-_tj[j]) < 1.0e-15 ) return 0.; //if the difference is really small

		den+=_betaj[j]/(t-_tj[j]);
	}
	
	double res=0.;
	res = _betaj[j]/(t-_tj[j]) / den;
	return res;
}
 
vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const
{
	vector_d fj(_p+1, 0.);;
	for(size_t i=0; i<=_p; i++){
		fj[i] = fnc(_tj[i]);
	}
	return fj;
}

double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const
{
	//domain checks
   if (t < -1 || t > 1 || (end - start) != _p) {
      throw std::domain_error("[StandardGrid::interpolate]: t=" + std::to_string(t)
                              + " \\notin [-1, +1] OR view "
                                "into fj of wrong size: ["
                              + std::to_string(start) + ", " + std::to_string(end) + "]");
   }

   	//compute the denominator
   double den = 0;
   for (size_t l = 0; l <= _p; l++) {
      if (fabs(t - _tj[l]) < 1.0e-15) {
         double sum = 0;
         for (size_t i = start, j = 0; i <= end; i++, j++) {
            sum += fj[i] * _Dij[j][l];
         }
         return sum;
      }
      den += _betaj[l] / (t - _tj[l]);
   }

   double sum = 0;
   for (size_t i = start, j = 0; i <= end; i++, j++) {
      sum += poli_weight_der(t, j, den) * fj[i];
   }
   return sum; //nelle note sarebbe p'(t)
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



double StandardGrid::poli_weight_der(double t, size_t j, double den) const //nelle note sarebbe lj'(t)
{
	//input checks
	if (t < -1 || t > 1) {
      throw std::domain_error("[StandardGrid::poli_weight]: t=" + std::to_string(t)
                              + " \\notin [-1, +1]");
   }

	double sum=0;
	for (size_t i = 0; i<= _p; i++){
		if(fabs(t-_tj[i]) < 1.0e-15) return _Dij[j][i]; //se valuto su un punto della griglia
		sum+=_Dij[j][i]*poli_weight(t, i, den);
	}
	return sum;
} 
 




void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const //nelle note: cambia fi con f_tildej
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






} // namespace Chebyshev
} // namespace Interpolation
