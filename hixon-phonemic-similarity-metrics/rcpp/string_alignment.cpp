#include <Rcpp.h>
#include <string.h>

//[[Rcpp::export]]
int string_dist(std::string s, std::string t) {
  int n = s.size();
  int m = t.size();
  Rcpp::IntegerMatrix d(n + 1, m + 1);
  
  if (n == 0) return n;
  
  if (m == 0) return m;
  
  for(int i = 0; i <= n; i++) {
    d(i, 0) = i;
  }
  
  for(int j = 0; j <= m; j++) {
    d(0, j) = j;
  }
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      if(s[i - 1] == t[j - 1]) {
        d(i, j) = d(i - 1, j - 1); // No operation! if the same character
      }
      else {
        d(i, j) = std::min(
          d(i - 1, j) + 1,
          std::min(
            d(i, j - 1) + 1,
            d(i - 1, j - 1) + 2
          )
        );
      }
    }
  }
  
  return d(n, m);
}

//[[Rcpp::export]]
Rcpp::IntegerMatrix string_matrix(std::string s, std::string t) {
  int n = s.size();
  int m = t.size();
  Rcpp::IntegerMatrix d(n + 1, m + 1);
  
  if (n == 0) return n;
  
  if (m == 0) return m;
  
  for(int i = 0; i <= n; i++) {
    d(i, 0) = i;
  }
  
  for(int j = 0; j <= m; j++) {
    d(0, j) = j;
  }
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      if(s[i - 1] == t[j - 1]) {
        d(i, j) = d(i - 1, j - 1); // No operation! if the same character
      }
      else {
        d(i, j) = std::min(
          d(i - 1, j) + 1,
          std::min(
            d(i, j - 1) + 1,
            d(i - 1, j - 1) + 2
          )
        );
      }
    }
  }
  
  return d;
}

//[[Rcpp::export]]
Rcpp::CharacterVector string_align(std::string s1, std::string s2) {
  
  int n = s1.size();
  int m = s2.size();
  
  Rcpp::IntegerMatrix d = string_matrix(s1, s2);
  
  std::string string_s1 ("");
  std::string string_s2 ("");
  
  int i = n;
  int j = m;
  int l = n + m;
  int s1pos = l;
  int s2pos = l;
  int temp = 0;
  int str_s1[l+1];
  int str_s2[l+1];
  
  while(i != 0 && j != 0) {
    if(s1[i - 1] == s2[j - 1]) {
      temp = d(i - 1, j - 1); // if same character then keep unchanged.
    }
    else temp = d(i - 1, j - 1) + 2; // assign substitution cost when different.
    
    if(d(i, j) == temp) {
      // move down diagonally when no substitution has been made.
      str_s1[s1pos--] = (int)s1[i - 1];
      str_s2[s2pos--] = (int)s2[j - 1];
      i--;
      j--;
    }
    else if(d(i, j) == 1 + d(i - 1, j)) {
      str_s1[s1pos--] = (int)s1[i - 1];
      str_s2[s2pos--] = (int)'*';
      i--;
    }
    else if(d(i, j) == 1 + d(i, j - 1)) {
      str_s1[s1pos--] = (int)'*';
      str_s2[s2pos--] = (int)s2[j - 1];
      j--;
    }
  }

  // Fill rest of it with *
  while(s1pos > 0) {
    if(i > 0) str_s1[s1pos--] = (int)s1[--i];
    else str_s1[s1pos--] = (int)'*';
  }
  
  while(s2pos > 0) {
    if(j > 0) str_s2[s2pos--] = (int)s2[--j];
    else str_s2[s2pos--] = (int)'*';
  }
  
  int ix = 1;
  
  // if both * then move the starting point ahead.
  for(int i = l; i > 0; i--) {
    if((char)str_s1[i] == '*' && (char)str_s2[i] == '*') {
      ix = i + 1;
      break;
    }
  }
  
  for(int i = ix; i <= l; i++) {
    string_s1 += (char)str_s1[i];
  }
  
  for(int j = ix; j <= l; j++) {
    string_s2 += (char)str_s2[j];
  }
  
  Rcpp::CharacterVector c = Rcpp::CharacterVector::create(string_s1, string_s2);
  
  return(c);
}
