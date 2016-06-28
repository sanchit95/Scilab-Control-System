// transfer function tf()
// SYS = TF(NUM,DEN) creates a continuous-time transfer function SYS with numerator NUM and denominator DEN.
// It works for SISO as well as for MIMO systems.
//SYS = TF(NUM,DEN,TS) creates a discrete-time transfer function with sampling time TS 
//(set TS = 'd' if the sampling time is undetermined).
//S = TF('s') specifies the transfer function variable.
//
//[sys details] = sys = tf([1 2],[1 0 10], 0.5) 
// Here details gives the structure containing all the system parameters such:-
//   SystemType: "rational"
//   num: [1 2]
//   den: [1 0 10]
//   Ts: 0.5
//   InputName: ""
//   OutputName: ""  
//
// Example 1
// ss = tf([3 2 3],[4 3 2 10],0.1) 
//
// Example 2
// a= tf({[3 42 -2];[-6 -84 -9]; [2 3 4]},{[3 4 5];[-3 -21 -7]; [5 3 4]})
//
// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/
// http://in.mathworks.com/help/control/ ; 
// https://en.wikipedia.org/wiki/Control_systems;
//
//
// Author (s):
// Sanchit Gupta
//-----------------------------------------------------------------------------------------------------------------------------//
function[sys, details] = tf(varargin)
    if (length(varargin)==2)&(typeof(varargin(1))=='constant')&(typeof(varargin(2))=='constant')
        // case of continuous time system
        s=poly(0,'s') ;
        num_full = varargin(1);
        den_full = varargin(2) ;
        if (size(varargin(1),1)<>size(varargin(2),1))
            error('Num and Den must be arrays of compatible sizes.')
        end
        for j = 1:size(num_full,1)
            num = num_full(j,:);
            den = den_full(j,:) ;
            if (den==0)|(num==[])|(den==[])
                error('Numerator must be a non-empty matrix. Denominator must be a non-empty and non-zero matrix.')
            end
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
            sys(j,1) = syslin('c',num2/den2) ;
        end 
        details = struct('SystemType',typeof(sys),'num',num_full,'den',den_full,'Ts',sys.dt,'InputName','','OutputName','') ;
    elseif (length(varargin)==3)&(typeof(varargin(1))=='constant')&(typeof(varargin(2))=='constant')& (typeof(varargin(3))=='constant')
        // case of discrete time system
        z = poly(0,'z') ;
        num_full = varargin(1);
        den_full = varargin(2) ;
        if (size(varargin(1),1)<>size(varargin(2),1))
            error('Num and Den must be arrays of compatible sizes.')
        end
        for j = 1:size(num_full,1)
            num = num_full(j,:);
            den = den_full(j,:) ;
            if (den==0)|(num==[])|(den==[])
                error('Numerator must be a non-empty matrix. Denominator must be a non-empty and non-zero matrix.')
            end
            ts = varargin(3) ;   // Sampling Time
            num2=0 ;
            den2=0 ;
            for i = 1:length(num)
                num1(i)=num(i)*z^(length(num)-i)
                num2=num2+num1(i)
            end
            for i = 1:length(den)
                den1(i)=den(i)*z^(length(den)-i)
                den2=den2+den1(i)
            end
            if(typeof(ts)=='string')&(ts=='d')
                sys(j,1) = syslin('d',num2/den2) ;
            elseif (size(ts)==1)&(ts>=0)
                sys(j,1) = syslin(ts,num2/den2) ;
            else
                error('Sampling time(ts) must be a positive scalar or equal to d to mean unspecified.') ;
            end
        end
        details = struct('SystemType',typeof(sys),'num',num_full,'den',den_full,'Ts',sys.dt,'InputName','','OutputName','') ;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    elseif (length(varargin)==3)&(typeof(varargin(1))=='constant')&(typeof(varargin(2))=='constant')& (typeof(varargin(3))=='string')&(varargin(3)=='c')
        // case of continuous time system
        s=poly(0,'s') ;
        num_full = varargin(1);
        den_full = varargin(2) ;
        if (size(varargin(1),1)<>size(varargin(2),1))
            error('Num and Den must be arrays of compatible sizes.')
        end
        for j = 1:size(num_full,1)
            num = num_full(j,:);
            den = den_full(j,:) ;
            if (den==0)|(num==[])|(den==[])
                error('Numerator must be a non-empty matrix. Denominator must be a non-empty and non-zero matrix.')
            end
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
            sys(j,1) = syslin('c',num2/den2) ;
        end 
        details = struct('SystemType',typeof(sys),'num',num_full,'den',den_full,'Ts',sys.dt,'InputName','','OutputName','') ;
    elseif (length(varargin)==3)&(typeof(varargin(1))=='constant')&(typeof(varargin(2))=='constant')& (typeof(varargin(3))=='string') & (varargin(3)=='d')
        // case of discrete time system
        z = poly(0,'z') ;
        num_full = varargin(1);
        den_full = varargin(2) ;
        if (size(varargin(1),1)<>size(varargin(2),1))
            error('Num and Den must be arrays of compatible sizes.')
        end
        for j = 1:size(num_full,1)
            num = num_full(j,:);
            den = den_full(j,:) ;
            if (den==0)|(num==[])|(den==[])
                error('Numerator must be a non-empty matrix. Denominator must be a non-empty and non-zero matrix.')
            end
            ts = varargin(3) ;   // Sampling Time
            num2=0 ;
            den2=0 ;
            for i = 1:length(num)
                num1(i)=num(i)*z^(length(num)-i)
                num2=num2+num1(i)
            end
            for i = 1:length(den)
                den1(i)=den(i)*z^(length(den)-i)
                den2=den2+den1(i)
            end
            if(typeof(ts)=='string')&(ts=='d')
                sys(j,1) = syslin('d',num2/den2) ;
            elseif (size(ts)==1)&(ts>=0)
                sys(j,1) = syslin(ts,num2/den2) ;
            else
                error('Sampling time(ts) must be a positive scalar or equal to d to mean unspecified.') ;
            end
        end
        details = struct('SystemType',typeof(sys),'num',num_full,'den',den_full,'Ts',sys.dt,'InputName','','OutputName','') ;
    elseif (typeof(varargin(1))=='string')&(length(varargin)==1) 
            sys = poly(0,varargin(1)) ;
    else
        error('Incorrect input arguments.') ;
    end
endfunction
