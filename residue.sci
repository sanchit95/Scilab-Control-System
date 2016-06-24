// residue function
//RESIDUE Partial-fraction expansion (residues).
//   [R,P,K] = RESIDUE(B,A) finds the residues, poles and direct term of
//   a partial fraction expansion of the ratio of two polynomials B(s)/A(s).
//   If there are no multiple roots,
//
//      B(s)       R(1)       R(2)             R(n)
//      ----  =  -------- + -------- + ... + -------- + K(s)    ;
//      A(s)     s - P(1)   s - P(2)         s - P(n)
//    
//[R,P] = RESIDUE(sys), where sys is in rational form finds 
// residue, poles of the system. 
//               B(s)
//      sys  =  ------  ;
//               A(s)  
//
// Example 1:
//[r p k] = residue([20],[2 3 4 9 8]) ;
//
// Example 2:
// s = %s ;
// sys = syslin('c',(s+4)/(s^2+5*s+10)) ;
//[r p] = residue(sys)
//
//Warning: Numerically, the partial fraction expansion of a ratio of
//polynomials represents an ill-posed problem.  If the denominator
//polynomial, A(s), is near a polynomial with multiple roots, then
//small changes in the data, including roundoff errors, can make
//arbitrarily large changes in the resulting poles and residues.
//Problem formulations making use of state-space or zero-pole
//representations are preferable.
//
//
// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/ ;
// http://in.mathworks.com/help/control/ ; 
// https://en.wikipedia.org/wiki/Control_systems ;
//
// Author (s):
// Sanchit Gupta
//
//-----------------------------------------------------------------------------------------------------------------------//

function[r, poles, k] = residue(varargin)
    temp = [] ;
    k = [] ;
    if (length(varargin)==1) & (typeof(varargin(1))=='rational')
        pf = pfss(varargin(1));
        temp = pf
        if typeof(pf(length(pf)))=='polynomial'
            k = pf(length(pf)) ;
        end
    elseif (length(varargin)==2) & (typeof(varargin(1))=='constant') & (typeof(varargin(2))=='constant')
        num = varargin(1) ;
        den = varargin(2) ;
        if den == 0
            error('The denominator must me nonzero.')
        end
        s=poly(0,'s') ;
        num2=0 ;
        den2=0 ;
        for i = 1:length(num)
            num1(i)=num(i)*s^(length(num)-i)
            num2=num2+num1(i)
        end
        for i = 1:length(den)
            den1(i)=den(i)*s^(length(den)-i)
            den2=den2+den1(i)
        end
        f = syslin('c',num2/den2) ; 
        pf = pfss(f);
        temp = pf
        if typeof(pf(length(pf)))=='polynomial'
            k = pf(length(pf)) ;
        end
    else
        error('Incorrect Input arguments') ;
    end
    //-------------------------------------------------------------------------------------------------------------------//
    r = [] ;
    poles = [] ;
    lenth = length(temp) ;
    if  (typeof(pf(length(pf)))=='polynomial')
        lenth = lenth - 1 ;
    end
    for j = 1:lenth ;
        //        if (typeof(temp(j).den)=='constant')& (typeof(coeff(temp(j).num))=='constant')
        //            k = coeff(temp(j).num)/coeff(temp(j).den) ;
        //            break ;
        //        end
        root = roots(temp(j).den) ;
        if length(root) <= 2
            if length(root)==1
                r = [r ;coeff(temp(j).num);] ;
                poles = [poles;root;] ;
            elseif length(root)==2 & (root(1)<>root(2))
                //
                // temp(k) in the form
                //     q + p*s
                // -------------
                // as^2 + b*s + c
                //
                // Where root(1) and root(2) gives the roots of the polyomial  
                //
                //Then its solution is given by 
                //       m                  n
                //  -----------    +   -----------  ;
                //   s -root(1)        s - root(2)
                //
                //Where 
                //          q + p*root(1)                   -q - p*root(2)
                // m = ---------------------    ;    n = ------------------- ;
                //       root(1) - root(2)                root(1) - root(2)
                //
                numerator = coeff(temp(j).num);
                denominator = coeff(temp(j).den) ;
                p =0 
                if (length(numerator)==2) ;
                    p = numerator(2);
                end
                q = numerator(1);
                a = denominator(3) ; b = denominator(2); c = denominator(3) ;
                m = (q+p*root(1))/(root(1)-root(2)) ;
                n = (-q-p*root(2))/(root(1)-root(2)) ;
                r = [r;m;n;] ;
                poles = [poles;root(1);root(2);] ;
            elseif length(root)==2 & (root(1) == root(2))
                r = [r ;0;coeff(temp(j).num);] ;
                poles = [poles;root(1);root(2)] ;
            end
            //--------------------------------------------------------------------------------------------------------------//
        else
            // length(root)> 2
            num = temp(j).num ;
            den = 1;
            for count = 1:length(root)
                // +(0.00001*count)*sign((-1)^(count))
                den = den*(s-root(count)+(0.00001*count)*sign((-1)^(count)));
            end
            pff = pfss(num/(real(den))) ;
            disp(pff) ;
            len = length(pff) ;
            if  (typeof(pf(length(pf)))=='polynomial')
                len = len - 1 ;
            end
            for pp = 1:len ;
                //        if (typeof(temp(j).den)=='constant')& (typeof(coeff(temp(j).num))=='constant')
                //            k = coeff(temp(j).num)/coeff(temp(j).den) ;
                //            break ;
                //        end
                root = roots(pff(pp).den) ;
                if length(root) <= 2
                    if length(root)==1
                        r = [r ;coeff(pff(pp).num);] ;
                        poles = [poles;root;] ;
                    elseif length(root)==2 & (root(1)<>root(2))
                        //
                        // temp(k) in the form
                        //     q + p*s
                        // -------------
                        // as^2 + b*s + c
                        //
                        // Where root(1) and root(2) gives the roots of the polyomial  
                        //
                        //Then its solution is given by 
                        //       m                  n
                        //  -----------    +   -----------  ;
                        //   s -root(1)        s - root(2)
                        //
                        //Where 
                        //          q + p*root(1)                   -q - p*root(2)
                        // m = ---------------------    ;    n = ------------------- ;
                        //       root(1) - root(2)                root(1) - root(2)
                        //
                        numerator = coeff(pff(pp).num);
                        denominator = coeff(pff(pp).den) ;
                        p =0 
                        if (length(numerator)==2) ;
                            p = numerator(2);
                        end
                        q = numerator(1);
                        a = denominator(3) ; b = denominator(2); c = denominator(3) ;
                        m = (q+p*root(1))/(root(1)-root(2)) ;
                        n = (-q-p*root(2))/(root(1)-root(2)) ;
                        r = [r;m;n;] ;
                        poles = [poles;root(1);root(2);] ;
                    elseif length(root)==2 & (root(1) == root(2))
                        r = [r ;0;coeff(temp(j).num);] ;
                        poles = [poles;root(1);root(2)] ;
                    end
                end 
            end
        end
    end
endfunction
