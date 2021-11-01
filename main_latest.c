#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

clock_t start_BM, finish_BM,start_KMP, finish_KMP, start_BF, finish_BF;


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MAXCHAR 211
int * text;
int * pattern;
int * BA; //Border Array
int * PA; //Prefix Array
int * indetlist; //list conisting of positions at which the pattern occurs in text
int * prime;
int *arr;
int alphabet_size [] = {2, 3, 5, 7, 6, 10, 14, 15, 21, 35, 30, 42, 70, 105, 210};
//int ** bad_char_table;
typedef  int (*bad_character_rule_ptr)(int i, int c, int* amap, int len_amap, int**  bad_char, int len_bad_char);
typedef  int (*good_suffix_rule_ptr)(int i, int* big_l_array, int len_big_l_array, int* small_l_prime_array, int len_small_l_prime_array);
typedef  int (*match_skip_ptr)(int* small_l_prime_array, int len_small_l_prime_array);

typedef  struct _boyer_moore {
    int*  amap_s;  /*  in python it is dict.  */
    int**  bad_char_s;
    int  len_amap_s;
    int  len_bad_char_s;
    int* big_l_prime_array_s;
    int* big_l_array_s;
    int* small_l_prime_array_s;
    bad_character_rule_ptr  bad_character_rule_ptr_fun;
    good_suffix_rule_ptr  good_suffix_rule_ptr_fun;
    match_skip_ptr  match_skip_ptr_fun;
    struct _boyer_moore* next_boyer_moore;
}boyer_moore;

int * readbin1(char* address);
void borderarray(int * pat, int n);
int gcdExtended(int a, int b, int *x, int *y);
void init_boyer_moore_struct(int* p, int len_p, boyer_moore* p_bm);

void gen_prime(int sigma);
int gcd(int a, int b);
void gen_prime(int sigma);
bool indet(int x, int sigma);
int min(int i, int j);
int match(int * text, int pos1, int pos2, int j, int n);
void prefixarray(int * text, int n);
int compute_shift(bool indettext, int i, int j, int m_ell);
void KMP_Indet(int n, int m, int sigma, int m_ell);
void bruteforce(int n, int m);
int* z_array(int* s, int len);
int* n_array(int* s, int len);
int big_l_prime_array(int* p, int* n, int len, int* big_l_prime_array_tmp);
int big_l_array(int* p, int* lp, int len, int* big_l_array_tmp);
int small_l_prime_array(int* n, int len, int* small_l_prime_array_tmp);
void good_suffix_table(int* p, int len, int* big_l_prime_array, int* big_l_array, int* small_l_prime_array);
int good_suffix_mismatch(int i, int* big_l_prime, int len, int* small_l_prime);
int good_suffix_match(int* small_l_prime, int len);
int dense_bad_char_tab(int* p, int len_p, int* amap, int len_amap, int** tab); /*  length of tab array is in function return  */
//void table_bad_char(int* pattern, int len_p);
//int bad_char_rule(int pos, int mismatch_char );
void calculate_table(int * pattern, int len_p, int alphabet_size);
//--------------------------
int bad_character_rule_fun(int i, int c, int* amap, int len_amap, int**  bad_char, int len_bad_char);
int indet_good_suffix(int* pattern, int len_p, int* text, int i, int j);
int good_suffix_rule_fun(int i, int* big_l_array, int len_big_l_array, int* small_l_prime_array, int len_small_l_prime_array);
int match_skip_fun(int* small_l_prime_array, int len_small_l_prime_array);
int boyer_moore_fun(int* p, int len_p, boyer_moore p_bm, int* t, int len_t, int* occurrences, int sigma); /*  length of occurrences array is in function return    */
//-------------------------

int * readbin1(char* address){
    int *A;
    FILE *fp;
    fp = fopen(address, "rb");
    fseek(fp,0L,SEEK_END);
    int n = ftell(fp);
    //printf("\n n: %d\n", n);
    n = n/4;
    rewind(fp);
    //printf("n=%d \n", n);
    A = (int *)malloc(sizeof(int)*n);
    fread(A, sizeof(int),n,fp);

    if(fp != NULL){
        fread(A, sizeof(int),n,fp);

    }

    fclose(fp);
    return A;
}


void gen_prime(int sigma){
    prime = (int *)malloc(sizeof(int)*sigma);
    int n = sigma;
    int temp=1;
    int i=3, count, c;
    if (n>=1){
        prime[0] = 2;
    }
    for(count=2; count<=n; i++){
        for(c=2; c<i; c++){
            if(i%c== 0)
                break;
        }

        if(c==i){
            prime[temp] = i;
            count++;
            temp++;
        }
    }
}

void borderarray(int * pat, int n){
    if(n==0){return;}
    int i, b;
    BA = (int *)malloc(sizeof(int)*n);
    BA[0] = 0;
    for(i=0; i < n-1; i++){
        b = BA[i];
        while(b> 0 && pat[i+1] != pat[b]){ //text[i+1], text[b] are compared as the indicies begin with 0 instead of 1.
            b = BA[b-1];
        }
        if(pat[b]==pat[i+1]){
            BA[i+1]=b+1;
        }
        else{
            BA[i+1] = 0;
        }
    }
}

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b%a, a);
}

int gcdExtended(int a, int b, int *x, int *y)
{
    // Base Case
    if (a == 0)
    {
        *x = 0;
        *y = 1;
        return b;
    }

    int x1, y1; // To store results of recursive call
    int gcd = gcdExtended(b%a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b/a) * x1;
    *y = x1;

    return gcd;
}

bool indet(int x, int sigma){ //function to test if a letter is indeterminate or not.
    int i;
    int n = sigma;
    for(i=0; i < n; i++){
        if(x==prime[i]){break;}
    }
    if(i<n){
        return false;}
    return true;
}


int match(int * text, int pos1, int pos2, int j, int n){
    int x, y; //variables for gcdExtended function.
    while(gcd(text[pos1], text[pos2])> 1 && pos1 < n && pos2 < n){
        pos1 = pos1 + 1;
        pos2 = pos2 + 1;
    }
    //printf("pos1=%d\n\n\n",pos1);
    return pos1;
}

int min(int i, int j){
    if(i > j){ i = j;}
    return i;
}

void prefixarray(int * text, int n){
    int x, y; //variables for gcdExtended function.
    int lambda = text[0];
    text[n] = -1;
    int i, j, lim, pref; //end;
    PA = (int *)malloc(sizeof(int)*(n));
    for(i=0;i<n;i++){
        if(gcdExtended(text[i],lambda, &x, &y) > 1){
            PA[i] = 1;
        }
        else{PA[i]=0;}// EQUAL[i]=false;}
    }
    PA[0]=n;
    j = 1; lim = 1;
    while(j < n){
        if(PA[j] != 0){
            pref = PA[j]; //1 - initially the current length, then the maximum length, of the substring beginning at x[i] that matches a prefix of x;
            pref = match(text, pref, j+pref, j, n); //match(1, 2, 0 )
            PA[j]=pref;
        }
        j = j +1;
    }
}

int compute_shift(bool indettext, int i, int j, int m_ell){
    int max = 0;
    int * patnew;
    patnew = (int *)malloc(sizeof(int)*(2*j));
    if(indettext || j > m_ell-1){
        memcpy(patnew, pattern, (j)*sizeof(int));
        memcpy(patnew+j, text+i-j+1, (j)*sizeof(int));
        prefixarray(patnew, 2*(j));
        for(i=j;i<=2*j;i++){

            if(2*j-i == PA[i] && max < PA[i]){
                max=PA[i];
            }
        }
        j = max;
        //printf("compute shift j=%d\n",j);
    }
    else{
        j = BA[j];
        //printf("compute shiftj=%d\n",j);

    }
    return j;
}

void KMP_Indet(int n, int m, int sigma, int m_ell){

    int x, y; //variables for gcdExtended function.
    int i, j;
    //int m_ell; //length of the longest regular prefix of the pattern.
    bool indettext = false;
    //compute m_ell
    /*
    m_ell=m; //assumes that pattern is regular
    i = 0;
    while(i<m){
        if(indet(pattern[i], sigma)){

            m_ell=i;
            break;
        }
        i=i+1;
    }
    */
    //borderarray(pattern, m_ell); //Border array of the longest regular prefix of pattern.
    i=-1;//reset
    j=-1;
    int count = -1;
    //i and j are the index positions in the text and pattern that match.
    //As a result the substrings text[i-j..i] and pattern[0..j] match.
    while(i < n-1){
        if(gcd(text[i+1], pattern[j+1]) > 1){
            //printf("text[%d] equals pattern[%d]\n", i+1, j+1);
            if(indet(text[i+1], sigma)){
                indettext=true;
            }
            j = j+1; //j is the index positions the pattern such that the prefix of length j+1 has match with a substsring of text.
            i = i+1;

            if(j==m-1){
                //printf("the pattern found in this position of text: %d\n", i-j);
                count = count + 1;

                indetlist[count] = i-j;
                //printf("count: %d\t", count);
                //printf("%d\n", indetlist[count]);
                j = compute_shift(indettext, i, j, m_ell)-1;
                //printf("shift=j=%d\n", j);
                indettext = false;
            }
        }
        else{
            if(j==-1){
                i = i+1;
            }
            else{
                j = compute_shift(indettext, i, j, m_ell)-1;
                //  printf("shift=j=%d\n",j);
                indettext = false;
            }
        }
    }
}


void bruteforce(int n, int m){
    //int x;
    //int y;//variables for gcdExtended function.
    //When the pattern is algined at 'index' position text[i-j..i] match with pattern[0..j]
    int i, j;
    //is the position in the text for which we are checking wheather the pattern matches or not when aligned at this position.
    int index = 0;
    j=-1;
    i= index-1;
    int count = -1;
    while(index <= n-m){
        if(gcd(text[i+1], pattern[j+1]) > 1){
            j = j+1;
            i = i+1;
            if(j==m-1){
                count = count + 1;
                indetlist[count] = index;
                //printf("count: %d\t", count);
                //printf("%d \n", indetlist[count]);
                index = index+1;
                i = index-1;
                j = -1;
            }
        }
        else{
            if(j==-1){
                index = index+1;
                i = index-1;
            }
            else{
                index = index+1;
                i = index-1;
                //i = i-j;
                j = -1;
            }
        }
    }

}

//Z(i) be the length of the longest substring of S that starts at i and matches a prefix of S.
int* z_array(int* s, int len) /*  s : array, l : length  */
{
    if(len <= 1){
        return NULL;
    }
    int* z;
    z = (int*) calloc(len, sizeof(int));
    z[0] = len;
    for(int i = 1; i < len; i++){
        if(gcd(s[i], s[i -1]) > 1){
            z[1] +=1;
        }
        else{
            break;
        }
    }
    int r = 0;
    int l = 0;
    if(z[1] > 0){
        r = z[1];
        l = 1;
    }
    for(int k = 2; k < len; k++){
        if(z[k] != 0){
            return NULL;
        }
        if(k > r){
            for(int i = k; i < len; i++){
                if(gcd(s[i], s[i - k]) > 1){
                    z[k] += 1;
                }
                else{
                    break;
                }
            }
            r = k + z[k] - 1;
            l = k;

        }
        else{
            int nbeta = r - k + 1;
            int zkp = z[k - 1];
            if(nbeta > zkp){
                z[k] = zkp;
            }
            else{
                int nmatch = 0;
                for(int i = r + 1; i < len; i++){
                    if(gcd(s[i], s[i - k]) > 1){
                        nmatch += 1;
                    }
                    else{
                        break;
                    }
                }
                l = k;
                r = r + nmatch;
                z[k] = r - k + 1;
            }

        }

    }
    return z;
}


int* n_array(int* s, int len) /*  s : array, l : length  */
{
    int *tmp_s;
    tmp_s = (int*)malloc(len * sizeof(int));
    for(int i = 0; i < len; i++){
        tmp_s[i] = s[len - 1 - i];
    }
    int *z_tmp = 0;
    int *z_tmp_res = 0;
    z_tmp = z_array(tmp_s,len);
    z_tmp_res = (int*)malloc(len * sizeof(int));
    for(int i = 0; i < len; i++){
        z_tmp_res[i] = z_tmp[len - 1 - i];
    }
    free(tmp_s);  /* should free memory  */
    free(z_tmp);
    return z_tmp_res;
}

int big_l_prime_array(int* p, int* n, int len, int* big_l_prime_array_tmp)
{
    int i = 0;
    for(int j = 0; j < len - 1; j++){
        i = len - n[j];
        if(i < len){
            big_l_prime_array_tmp[i] = j + 1;
        }
    }
    return 0;

}


// big_l_array_tmp = l array in Python code
int big_l_array(int* p, int* lp, int len, int* big_l_array_tmp)
{
    big_l_array_tmp[1] = lp[1];
    for(int j = 2; j < len; j++){
        big_l_array_tmp[j] = MAX(big_l_array_tmp[j - 1], lp[j]);
    }
    return 0;
}

// small_l_prime_array_tmp is like small_lp in python code
int small_l_prime_array(int* n, int len, int* small_l_prime_array_tmp)
{
    for(int j = 0; j < len; j++){
        if(n[j] == j + 1){
            small_l_prime_array_tmp[len - j - 1] = j + 1;
        }
    }
    for(int j = len - 2; j > -1; j--){
        if(small_l_prime_array_tmp[j] == 0){
            small_l_prime_array_tmp[j] = small_l_prime_array_tmp[j + 1];
        }
    }
    return 0;
}

void good_suffix_table(int* p, int len, int* big_l_prime_array_ptr,
                       int* big_l_array_ptr, int* small_l_prime_array_ptr) /* *_ptrs are outputs  */
                       {
    int *n = 0;
    n = n_array(p, len);
    int res1 = big_l_prime_array(p, n, len, big_l_prime_array_ptr);
    int res2 = big_l_array(p, big_l_prime_array_ptr, len, big_l_array_ptr);
    int res3 = small_l_prime_array(n, len, small_l_prime_array_ptr);

                       }
                       // len = length of small_l_prime in python code
    int good_suffix_mismatch(int i, int* big_l_prime, int len, int* small_l_prime)
                       {
    if(i > len){
        return -1;
    }
    if(i == len - 1){
        return 0;
    }
    i += 1;
    if (big_l_prime[i] > 0){
        return len - big_l_prime[i];
    }
    return len - small_l_prime[i];

                       }

                       // len = length of small_l_prime in python code
                       int good_suffix_match(int* small_l_prime, int len)
                       {
    return len - small_l_prime[1];
                       }

int dense_bad_char_tab(int* p, int len_p, int* amap, int len_amap, int** tab)
                       {
    int* nxt;
    int c = 0;
    int tab_len = 0;
    nxt = (int*)calloc(len_amap, sizeof(int));
    //printf("IN dense bad char FUN\n");
    for(int i = 0; i < len_p; i++){
        c = p[i];
        int res_in_amap = 0; /*  to see if c is in amap array  */
        for(int j = 0; j < len_amap; j++){
            if(c ==  amap[j]){
                res_in_amap = 1;
            }
        }
        if(res_in_amap == 0){
            //printf("c is not in amap\n");
            return -1;
        }
        tab_len++;
        for(int j = 0; j < len_amap; j++){
            tab[tab_len - 1][j] = nxt[j];
        }

        nxt[amap[c]] = i + 1;
    }
    free(nxt);
    return tab_len;
                       }
int bad_character_rule_fun(int i, int c, int* amap, int len_amap, int**  bad_char, int len_bad_char){
    //printf("IN bad_character_rule_fun FUN\n");
    int res_in_amap = 0; /*  to see if c is in amap array  */

    for(int j = 0; j < len_amap; j++){
        if(c == amap[j]){
            res_in_amap = 1;
        }
    }

    if(res_in_amap == 0){
        //printf("IN bad_character_rule_fun FUN result is not in amap c = %d \n", c);
        return -1;
    }
    if(i >= len_bad_char){
        //printf("error in bad_character_rule_fun");
        return -1;
    }
    int ci = amap[c];
    return i - (bad_char[i][ci] - 1);

}

void calculate_table(int * pattern, int len_p, int alphabet_size){

    //int len_alphabet = 15;
    arr = (int *)malloc(MAXCHAR * len_p * sizeof(int));
    int i, j;
    int row = MAXCHAR;
    int col = len_p;

    for (i = 0; i < row; i++){
        for (j = 0; j < col; j++){
            *(arr + i*col + j) = -1;
        }
    }

    for (i = 0; i < MAXCHAR; i++) {
        for (j = 0; j < col; j++) {
            if(gcd(i, pattern[j]) > 1){
                *(arr + i*col + j) = j;
            }
        }
    }
    
    for (i = 0; i < MAXCHAR; i++) {
        for (j = 1 ; j < col; j++) {
            if(*(arr + i*col + j - 1) > *(arr + i*col + j)){
                *(arr + i*col + j) = *(arr + i*col + j - 1);
            }
        }

    }
    /*
    printf("TABLE");
    for (i = 0; i < row; i++){
        for (j = 0; j < col; j++){
            printf("%d, ",*(arr + i*col + j));
        }
        printf("\n");
    }
    printf("END_TABLE");
    */
}

int good_suffix_rule_fun(int i, int* big_l_array, int len_big_l_array,
                         int* small_l_prime_array, int len_small_l_prime_array)
                         {
    //printf("IN good_suffix_rule_fun FUN len_big_l_array = %d\n", len_big_l_array);
    if(i >= len_big_l_array){
        //printf("i > length of big_l_array\n");
        return -1;
    }
    if(i == len_big_l_array - 1){
        return 0;
    }
    int j = i + 1;
    if(big_l_array[j] > 0){
        return len_big_l_array - big_l_array[j];
    }
    return len_big_l_array - small_l_prime_array[j];

                         }

 int indet_good_suffix(int* pattern, int len_p, int* text, int i, int j){
    int* q_prime;
     int t_prime; //length of the matched substring to consider in the text
     //We consider the whole match substring in the text if its length is less than the pattern.
     //Otherwise, we consider its suffix of length j-1.
     if(j==len_p){t_prime=j-1;}
     else{t_prime=j;}
     
     //q_prime_length is the length of the q'. Its prefix is the reverse of match substring of length t_prime, and its suffix is the reverse pattern of length m-1.
     //We consider the suffix of length m-1 as we want to make sure that any alignments are not missed by only considering the prefix of pattern of length m-j as in the case of classical BM.
     int q_prime_length = t_prime + len_p - 1;
     q_prime = (int*)calloc(q_prime_length, sizeof(int));
     
     for(int l = 0; l < t_prime; l++){
         q_prime[l] = text[i + len_p-1 - l];
         //printf("q_prime[%d] = %d\n", l, text[i + len_p - 1 - l]);
     }
     //printf("\n");
     for(int l = 0; l < len_p-1; l++){
          q_prime[t_prime + l] = pattern[len_p-2-l];
         //printf("q_prime[%d] = %d\n", t_prime +l, pattern[len_p-2-l]);
      }
    
    //Construct the prefix array of q'
    prefixarray(q_prime, q_prime_length);
     for(int l = 0; l < q_prime_length; l++){
         //printf("PA[%d] = %d ", l, PA[l]);
      }
     //printf("\n");

     int gs_shift = 0; //It is the shift computed by adopting the good suffix rule.
     int rindex = t_prime; //It is the index position of the right most occurrence/match of t' (or its longest suffix) in p[1..m-1]
     int k = t_prime+1;
     while (k < q_prime_length){ //we traverse the prefix array to compute r-index from t-prime index to the last one, to ensure that
         //We only change the r-index if the prefix array value at the current index k is greater, as we are in search of longest suffix/match of t' and the rightmost occurs of the same.
         if (PA[k] > PA[rindex]){
             rindex = k;
             //printf("rindex: %d\n", rindex);
             
             //The below if condition ensures that the right most occurrence of t' is considered. We return the gs_shift here as we are not interested in a match > |t'| as it would include the characters from the pattern.
             if(PA[rindex] == t_prime){
                 //gs_shift returns the shift or the numbers of alignments to skip so that t'
                 gs_shift = len_p + rindex - q_prime_length;
                 //printf("gs_shift 1: %d, i=%d, len_p=%d, rindex=%d, q_prime_length = %d\n", gs_shift, i , len_p,  rindex, q_prime_length);
                 return gs_shift;
             }
              
         }
         k++;
     }
     //printf("q_prime_length: %d, rindex: %d\n", q_prime_length, rindex);
     //If no match is found, we shift the pattern past t' in the text, else if the match found has length < |t'| we return the the shift or the numbers of alignments to skip so that the rightmost longest matched suffix is appropriately aligned under t' in the text.
     if(PA[rindex] ==0){
         gs_shift = t_prime+1;
     }
     else
     {//gs_shift = q_prime_length-rindex-PA[rindex];
         gs_shift = len_p + rindex - q_prime_length;
     }
     //printf("gs_shift 2: %d\n", gs_shift);
     return gs_shift;
}

int match_skip_fun(int* small_l_prime_array, int len_small_l_prime_array)
{
    return len_small_l_prime_array - small_l_prime_array[1];
}
/*
void table_bad_char(int* pattern, int len_p){
    int bad_char_table [MAXCHAR][len_p];
    for (int i = 0; i < MAXCHAR; i ++){
        for(int j = 0; j < len_p; j++){
            if (gcd(pattern[i], pattern[j]) > 1){
                bad_char_table[i][j] = j;
            }
            else{
                bad_char_table[i][j] = bad_char_table[i][j-1];
            }
        }

    }


}

int bad_char_rule(int pos, int mismatch_char ){
    return pos - bad_char_table[mismatch_char][pos];
}
 */

int boyer_moore_fun(int* p, int len_p, boyer_moore p_bm, int* t, int len_t, int* occurrences, int sigma)
{
    int* p_reverse;
    p_reverse = (int*)calloc(len_p, sizeof(int));
    for (int k = 0; k < len_p; k++){
        p_reverse[k] = p[len_p - k - 1];
    }

    int i = 0;
    int shift = 0;
    int m_ell_suf = 0;
    bool indet_letter = false;
    bool mismatched = false;

    int l = 0;
    //P should be reverse
    //length of the longest regular prefix of P reverse = l
    while(l < len_p){
        if(indet(p_reverse[l], sigma)){
            m_ell_suf=l;
            //printf("m_ell-suf = %d", m_ell_suf);
            break;
        }
        l = l+1;
    }
    //printf("the m_ell_suf for %d = %d", i, m_ell_suf);
    if( l == len_p){
        m_ell_suf = l;
    }

    int occurrences_len = 0;
    int skip_bc = 0;
    int skip_gs = 0;
    while(i < len_t - len_p + 1){
        shift = 1;
        mismatched = false;
        for(int j = len_p - 1; j > -1; j--){
            if(indet(t[i+j], sigma)){
                indet_letter = true;
            }
            if(gcd(p[j], t[i + j]) == 1){
                //skip_bc = bad_character_rule_fun(j, t[i + j], p_bm.amap_s, p_bm.len_amap_s, p_bm.bad_char_s, p_bm.len_bad_char_s);
                skip_bc =  *(arr + t[i + j]*MAXCHAR + j);
                //printf("i=%d, j=%d, skip_bc = %d\n", i, j, skip_bc);
                skip_gs = 0;
                if(indet(t[i + j], sigma) || len_p - j - 1 >= m_ell_suf){
                    skip_gs = indet_good_suffix(p, len_p, t, i, len_p - j - 1);
                }
                else{
                    skip_gs = good_suffix_rule_fun(j, p_bm.big_l_array_s, len_p, p_bm.small_l_prime_array_s, len_p);
                }
                //printf("skip_gs=%d\n", skip_gs);
                shift = MAX(shift, skip_bc);
                shift = MAX(shift, skip_gs);
                mismatched = true;
                break;
            }

        }
        if(!mismatched){
            //printf("the index of match = %d\n", i);
            occurrences_len++;
            occurrences[occurrences_len - 1] = i;
            //printf("the index of match = %d\n", occurrences[occurrences_len - 1]);
            int skip_gs = 0;
            if(indet(t[i + len_p], sigma) ||  m_ell_suf < len_p){
                skip_gs = indet_good_suffix(p, len_p, t, i, len_p);
            }
            else{
                skip_gs = match_skip_fun(p_bm.small_l_prime_array_s, len_p);
            }

            shift = MAX(shift, skip_gs);
        }
       //printf("i=%d, shift = %d\n", i, shift);
        i += shift;
    }
    /*
    for(int i = 0; i < occurrences_len; i++){
        printf("[%d] : %d ", i, occurrences[i]);
    }
    */
    return occurrences_len;
}

void init_boyer_moore_struct(int* p, int len_p, boyer_moore* p_bm)
{


    p_bm->len_amap_s = 26;
    p_bm->amap_s = (int*)calloc(p_bm->len_amap_s, sizeof(int));  /* allocate nmap 0-25  */
    p_bm->big_l_prime_array_s = (int*)calloc(len_p, sizeof(int));
    p_bm->big_l_array_s = (int*)calloc(len_p, sizeof(int));
    p_bm->small_l_prime_array_s = (int*)calloc(len_p, sizeof(int));
    for(int j = 0; j < p_bm->len_amap_s; j++){
        p_bm->amap_s[j] = j;
    }
    p_bm->bad_char_s = (int**)malloc(len_p * sizeof(int*));
    for (int i = 0; i < len_p; i++){
        p_bm->bad_char_s[i] = (int*)malloc(p_bm->len_amap_s * sizeof(int));
    }
    //printf("IN init FUN AFTER amap_s allocation len amap = %d\n", p_bm->len_amap_s);
    p_bm->len_bad_char_s = dense_bad_char_tab(p, len_p, p_bm->amap_s, p_bm->len_amap_s, p_bm->bad_char_s);
    good_suffix_table(p, len_p, p_bm->big_l_prime_array_s, p_bm->big_l_array_s, p_bm->small_l_prime_array_s);

    p_bm->next_boyer_moore = NULL;

}
int main()
{
    int temp=0,n=0,m=0;
    char filename[100];
    char filename2[100];
    int *value;
    char* token;
    long int len;

    // read text file and generate .bin file of text
    int sigma = 4;//atoi(argv[1]);
    strcpy(filename, "demofile_text");
    FILE *text_file = fopen(filename, "r");
    fseek(text_file, 0, SEEK_END);
    len = ftell(text_file);
    rewind(text_file);
    char line1[len];
    strcpy(filename2, "demofile_text");
    strcat(filename2,"c.bin");
    FILE *text_fp = fopen(filename2, "wb");
    fread(line1, sizeof(char), len, text_fp);
    while (fgets(line1, len, text_file)){
        token = strtok(line1, ",");
        while (token != NULL) {
            n = n + 1;

            temp = atoi(token);
            value = &temp;
            fwrite (value, sizeof(int),1,text_fp);
            token = strtok(NULL, ",");
        }
    }

    fclose(text_file);
    fclose(text_fp);
    text = readbin1(filename2);

    strcpy(filename, "demofile_pattern");
    FILE *pattern_file = fopen(filename, "r");
    fseek(pattern_file, 0, SEEK_END);
    len = ftell(pattern_file);
    rewind(pattern_file);
    char line[len];
    strcpy(filename2, "demofile_pattern");
    strcat(filename2,"c.bin");
    FILE *pattern_fp = fopen(filename2, "wb");
    fread(line, sizeof(char), len, pattern_fp);
    while (fgets(line, len, pattern_file)){
        token = strtok(line, ",");
        while (token != NULL) {
            m = m + 1;

            temp = atoi(token);
            value = &temp;
            //printf("%d ", temp);
            fwrite (value, sizeof(int),1,pattern_fp);
            //printf("\n%d,", temp);
            token = strtok(NULL, ",");
        }
    }
    fclose(pattern_file);
    fclose(pattern_fp);
    pattern = readbin1(filename2);
    gen_prime(sigma);

    if (m>n) {
        fprintf(stderr,"Length of the pattern cannot be greater than text length\n");
        return -1;
    }
    indetlist = (int *)malloc(sizeof(int)*(n-m+1));
    //printf("-------Brute Force-------\n");
    start_BF = clock();
    bruteforce(n,m);
    finish_BF = clock();
    double BF_time = (double)((finish_BF - start_BF) / (double)CLOCKS_PER_SEC);
    printf("%.7f \n", BF_time);
    //printf("-------KMP Indet-------\n");
    int m_ell=m; //assumes that pattern is regular
    int i = 0;
    while(i<m){
        if(indet(pattern[i], sigma)){

            m_ell=i;
            break;
        }
        i=i+1;
    }

    borderarray(pattern, m_ell);
    start_KMP = clock();
    KMP_Indet(n,m, 4, m_ell);
    finish_KMP = clock();
    double KMP_time = (double)((finish_KMP - start_KMP) / (double)CLOCKS_PER_SEC);
    printf("%.7f \n", KMP_time);

    //int p[6] = {15, 6, 6, 3, 30, 6};
    //int t[24] = {7, 2, 15, 7, 5, 3, 7, 10, 30, 105, 6, 6, 3, 15, 6, 30, 7, 105, 15, 6, 3, 3, 2, 6};
    //printf("-------Boyer Moore-------\n");
    calculate_table(pattern, m, MAXCHAR);
    boyer_moore my_boyer_moore;
    init_boyer_moore_struct(pattern, m, &my_boyer_moore);
    int* occurrences = (int*)calloc(n, sizeof(int));
    start_BM = clock();
    int res = boyer_moore_fun(pattern, m, my_boyer_moore, text, n, occurrences, sigma);
    finish_BM = clock();
    double BM_time = (double)((finish_BM - start_BM) / (double)CLOCKS_PER_SEC);
    printf("%.7f \n", BM_time);
    //printf("IT is runing\n");
    //printf("res is: %d", res);
    /*
    for(int i = 0; i < n; i++){
        printf("%d ", text[i]);
    }
    printf("\n");
    for(int i = 0; i < m; i++){
        printf("%d ", pattern[i]);
    }
    */
    free(arr);
    //printf("----------------------------------------");
    return 0;
}





