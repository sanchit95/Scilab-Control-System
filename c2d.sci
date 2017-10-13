//SYSD = c2d(SYSC,TS,'METHOD')
//computes the discrete time model of the continuous time system SYSC, with the sampling time TS.
//The string METHOD selects the discretization method amongst the folowing:
//    'zoh'       Zero-order hold on the inputs 
//    'foh'       Linear interpolation of inputs
//    'impulse'   Impulse-invariant discretization
//    'tustin'    Bilinear (Tustin) approximation
//    'matched'   Matched pole-zero method (for SISO systems only)
//The default is 'zoh' when METHOD is omitted. The sampling time TS should 
//be specified in the time units of SYSC. 
//
//For setting prewarping frequency for Bilinear Tustin Method or Tustin method, use the function in the following way
//c2d(sysc,t,'tustin','PrewarpFrequency',.5);
//
// Example:
// s = %s
// sys = syslin('c',(s+5),(s^2+6*s+8));
// sysd = c2d(sys,0.2,'matched')
//
// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/ ;
// http://in.mathworks.com/help/control/ ; 
// https://en.wikipedia.org/wiki/Control_systems ;
//
// Author (s):
// Sanchit Gupta
//-----------------------------------------------------------------------------------------------------------------------------//


//Bilinear Tustin method
function [s1]=tustin(s,t,fp)
    //Syntax : s1=cls2dls(sl,t [,fp])
    //
    // Given sl=[a,b,c,d] (syslin list) continuous time system, cls2dls
    // returns the sampled system obatined by the bilinear transform
    // s=(2/t)*(z-1)/(z+1).
    //
    // t sampling period
    // fp prevarping frequency in hertz
    //!

    [lhs,rhs]=argn(0)
    if typeof(s)<>"state-space"
        s = tf2ss(s) ;
    end
    fs=1/t
    if rhs==3 then fp=2*%pi*fp;fs=fp/tan(fp/fs/2)/2,end //prewarping

    a=2*fs*eye()-s(2)
    ad=a\(2*fs*eye()+s(2))
    b=(ad+eye())/a*s(3);
    d=s(5)+s(4)/a*s(3)
    s1=ss2tf(syslin("d",ad,b,s(4),d,s(6))) ;
endfunction

//-----------------------------------------------------------------------------------------------------------------------------//
//Zero Order Hold
function[sysd] = zoh(sysc,t);
    s = %s ;
    z = %z ;
    if typeof(sysc)=="state-space"
        sysc = ss2tf(sysc) ;
    end
    [m n] = size(sysc) ;
    for i=1:m
        for j = 1:n
            z = roots(sysc(i,j).num);
            p = roots(sysc(i,j).den) ;
            if (length(z)>length(p))
                error('The zoh, foh, and impulse methods cannot be used for models with more zeros than poles.') ;    
            end
            sysd(i,j) = ss2tf(dscr(sysc(i,j),t)) ;
        end
    end
endfunction

//---------------------------------------------------------------------------------------------------------------------------//

//First Order Hold
function[sysd] = foh(sysc,t);
    s = %s ;
    z = %z ;
    if typeof(sysc)=="state-space"
        sysc = ss2tf(sysc) ;
    end
    sys_temp =(sysc)*(syslin('c',(t*s+1)/(t*s)));
    sys_temp1 = ss2tf(dscr(sys_temp,t)) ;
    sysd = (z-1)*sys_temp1 ;
    sysd = syslin(t,sysd) ;
endfunction

//----------------------------------------------------------------------------------------------------------------------------//

//Impulse invariant discretization
// If form of function is   b/(s+a) ;
//Then its z transform is bT/(1 - (e^(-aT))*z^-1)
// This method is only possible for a artictly proper transfer function(more poles than zeros).
function[sysd] = impulse_samp(sysc,t);
    s = %s ;
    z = %z ;
    if typeof(sysc)=="state-space"
        sysc = ss2tf(sysc) ;
    end
    [m n] = size(sysc)
    for i=1:m
        for j = 1:n
            z = roots(sysc(i,j).num);
            p = roots(sysc(i,j).den) ;
            if (length(z)>length(p))
                error('The zoh, foh, and impulse methods cannot be used for models with more zeros than poles.') ;    
            end
            temp = pfss(sysc(i,j));
            tempd = [] ;
            for k = 1:length(temp) 
                root = roots(temp(k).den) ;
                if length(root)==1
                    z = %z;
                    num = coeff(temp(k).num) ;
                    tempd(1,k) = syslin(t,(t*(num))/(1-exp(root*t)*(z^-1))) ;
                elseif length(root)==2 
                    //
                    // temp(k) in the form
                    //     p*s + q
                    // -------------
                    // as^2 + b*s + c
                    //
                    //Then its solution is given by 
                    //                  m                                      n
                    //  --------------------------------    +   --------------------------------  ;
                    //  s + b/(2*a) + sqrt(b^2-4ac)/(2a)        s + b/(2*a) - sqrt(b^2-4ac)/(2a)
                    //
                    //Where 
                    //      p*sqrt(b^2-4ac) + pb - 2aq            p*sqrt(b^2-4ac) - pb + 2aq
                    // m = ----------------------------   ;  n = ---------------------------- ;
                    //          2*sqrt(b^2-4ac)                        sqrt(b^2-4ac)
                    // 
                    if root(1)==root(2) 
                        root(2) = root(2)+0.00001 ;
                    end
                    z = %z ;
                    numerator = coeff(temp(k).num);
                    denominator = coeff(temp(k).den) ;
                    p =0 
                    if (length(numerator)==2) ;
                    p = numerator(2);
                    end
                    q = numerator(1);
                    a = denominator(3) ; b = denominator(2); c = denominator(3) ;
                    m = (q+p*root(1))/(root(1)-root(2)) ;
                    n = (-q-p*root(2))/(root(1)-root(2)) ;
                    sumd = syslin(t,(t*m)/(1-exp((root(1))*t)*(z^-1))) + syslin(t,(t*n)/(1-exp((root(2))*t)*(z^-1))) ;
                    tempd(1,k) = sumd 
                end                 
            end
            sysd(i,j) =sum(tempd);
        end
    end
    sysd = syslin(t,sysd) ;
endfunction

//-----------------------------------------------------------------------------------------------------------------------------//

//Matched pole zero method
//Poles and zeros of the continuous function will also be the poles and zeros of the discrete function.
// Put Z = exp(s*t) ;
//After mapping alll poles and zeros in the z domain, then do gain matching bu puttting
//         lt s->0 G(s)
//   K =  --------------
//         lt z->1 G(z)
//
function[sysd] = matched(sysc,t);
    s = %s ;
    z = %z ;
    if typeof(sysc)=="state-space"
        sysc = ss2tf(sysc) ;
    end
    [m n] = size(sysc)
    for i=1:m
        for j = 1:n
            zeroes = roots(sysc(i,j).num);
            poles = roots(sysc(i,j).den) ;
            
            if zeroes == []
                temp_num = 1;
            else
                temp_num = 1 ;
                for k = 1:length(zeroes) 
                    temp_num = temp_num*(z-exp(t*zeroes(k))) ;
                end    
            end
            if poles == []
                temp_den = (1+z) ;
            else
                temp_den = 1 ;
                for k = 1:length(poles) 
                temp_den = temp_den*(z-exp(t*poles(k))) ;
                end
            end
            temp = syslin(t,temp_num/temp_den) ;
            gain = horner(sysc(i,j),0.00001)/horner(temp,1.00001) ; 
            sysd(i,j) = gain*temp ;
        end
    end 
endfunction

//----------------------------------------------------------------------------------------------------------------------------//

// Main Function
function[sysd] = c2d(varargin)
    if length(varargin) < 2
        error('Not enough input arguments') ;
    end
    sysc = varargin(1) ;
    t = varargin(2) ;
    if (typeof(t)<> 'constant')|(length(t) <> 1) | (t <= 0) | (t == %inf)
        error('Time must be a scalar and positive real quantity.') ;
    end
    if (typeof(sysc)<>"state-space") & (typeof(sysc)<>"rational")
        error('Linear State space or Rational system expected') ;
    end
    if sysc.dt==[] then
        warning("Input argument is assumed continuous time.");
        sysc.dt="c";
    end
    if sysc.dt<>"c" then
        error("Wrong value for input argument. Continuous time system expected.") ;
    end
    if (length(varargin)==2)
        sysd = zoh(sysc,t);
    elseif (length(varargin)==3) & (typeof(varargin(3))=='string')
        method = varargin(3);
        if method == 'zoh'
            sysd = zoh(sysc,t);
        elseif method == 'foh'
            sysd = foh(sysc,t) ;
        elseif method == 'impulse'
            sysd = impulse_samp(sysc,t) ;
        elseif method == 'tustin'
            sysd = tustin(sysc,t) ;
        elseif method == 'matched'
            sysd = matched(sysc,t) ;
        else
            error("The discretization method of the c2d command must be one of the following strings: zoh, foh, impulse, tustin, or matched.");
        end
    elseif (length(varargin)==5)&(varargin(3)=='tustin')&(varargin(4)=='PrewarpFrequency')&(typeof(varargin(5))=='constant')& (size(varargin(5))==[1 1]) 
        fp = varargin(5) ;
        sysd = tustin(sysc,t,fp) ;
    else
        error('Incorrect input arguments') ;
    end
endfunction
    

