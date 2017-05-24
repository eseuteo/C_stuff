
#define LENGTH 16

int main(){
  int vector1[LENGTH] = {29, 44, 107, 37, 98, 175, 30, 106, 73, 114, 134, 51, 17, 76, 101, 38};
  int vector2[LENGTH];
  int k = LENGTH;
  double res;
  res = euclidean_distance(vector1, vector2, k);
}

double euclidean_distance(int * vector1, int * vector2, int k){
  double res;
  for (int i=0; i<k; i++){
    res += pow((vector1[i]-vector2[i]), 2);
  }
  return sqrt(res);
}
