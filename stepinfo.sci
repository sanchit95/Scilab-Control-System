//function stepinfo(sys)
//This function computes sytem characteristics when step input is given to the system. 
//stepinfo(Y,T,YFINAL) computes the following parameters assuming zero initial conditions and no input offset.
//  RiseTime: rise time of the system
//  SettlingTime: settling time
//  SettlingMin: min value of Y once the response has risen
//  SettlingMax: max value of Y once the response has risen
//  Overshoot: percentage overshoot (relative to YFINAL)
//  Undershoot: percentage undershoot
//  Peak: peak absolute value of Y
//  PeakTime: time at which peak absolute value is reached.
//
//
//stepinfo(Y,T,YFINAL)
//For SISO systems Y and T are vectors of same length NS. 
//For MIMO systems with NU inputs, NY outputs, Y should be 
//specified as an array of size NS-by-NY-by-NU. 
//Similarly, YFINAL should be specified as a NY-by-NU array. 
//Then the function returns an output array of each quantity 
//of size NY-by-NU.
//
//stepinfo(Y,T) uses the last sampled value of Y as the YFINAL
//stepinfo(Y) assumes T = 1:NS.
//
//stepinfo(sys) calculates the time response characteristics of 
//dynamic system sys. The Rise Time, Settling Time etc. are all 
//calculated in Time Units of the system.
//
//
//The default value of RiseTimeLims = [0.1 0.9]  ie. from 10% to 90% of the final value.
//The user may change the limits by setting RiseTimeLims to the desired values.
//Example :
//If the desired value of rise time is from 5% to 95% of the final value, then the command looks like 
//stepinfo(y,t,yfinal,'RiseTimeLims'[0.05 0.95])
//
//Settling time with default 2% limits
//The user may change this by setting SettlingTimeThreshold to the desired value.
//Example :
//If the desired value of settling time threshold is 5% of the final value, then the command looks like 
//stepinfo(sys,'SettlingTimeThreshold',0.05)  
//
//The function is also capable of giving Rise time high and rise time low. These values are not printed by default but they
//can only be assigned to a variable. 
//RiseTimeHigh, RiseTimeLow can be extracted using the function in the following way
//Example :
//[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10] = stepinfo(Y,T,YFINAL)
//Then the variables x9 and x10 will give the values of RiseTimeHigh and RiseTimeLow respectively.
//
//Example:
// s = %s ;
// sys = syslin('c',(s+6)/(s^2+6*s+4)) ;
// stepinfo(sys) ;
//
//
// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/
// http://in.mathworks.com/help/control/ ; 
// http://in.mathworks.com/help/control/ref/stepinfo.html ;
// https://en.wikipedia.org/wiki/Control_systems;
// 
//
// Author(s):
// Sanchit Gupta
//----------------------------------------------------------------------------------------------------------------------------//
//
//
//stepinfo(y,t,yfinal) 
//Peak Value, Peak Time 
function[peak_value, peak_time] = peak(y,t,yfinal) ;
    if abs(yfinal) < %inf
        [peak_value, index] = max(y) ; 
        peak_time = t(index) 
    else
        // System is not stable
        peak_value = %inf
        peak_time = %inf 
    end
endfunction

//Function to find maximum and minimum overshoot
function[overshoot, undershoot] = overshoot_and_undershoot(y,t,yfinal)
    if abs(yfinal) < %inf
        if yfinal==0
            overshoot = %inf;
            if y(1:length(y))>= 0
                undershoot = 0;
            else
                undershoot = %inf;
            end
       else
          yrel = y/yfinal;
          overshoot = 100 * max(0,max(yrel-1));
          undershoot = -100 * min(0,min(yrel));
       end 
    else
        // System is not stable
        overshoot = %nan
        undershoot = %nan
    end
endfunction


// Function to find Rise Time, settling minimum, settling maximum
// The default value of RiseTimeLims = [0.1 0.9]  ie. from 10% to 90% of the final value.
// The user may change the limits by setting RiseTimeLims to the desired values.
// Settling min and Settling maximum are the minimum and maximum values once the response has risen.
function[rise_time, settling_min, settling_max, tHigh,tLow] = risetime(y,t,yfinal,RiseTimeLims)
    if abs(yfinal) < %inf
        ns = length(t) ;
        yLow = y(1) + RiseTimeLims(1)*(yfinal - y(1)); 
        iLow = 1 + find((y(1:ns-1)-yLow).*(y(2:ns)- yLow) <= 0,1);
        if isempty(iLow) 
            // Has not reached the lower 
            tLow = %nan
        elseif (iLow > 1) & (y(iLow)~=y(iLow-1))
            // Interpolate for more accuracy
            tLow = t(iLow) + (t(iLow)-t(iLow-1))/(y(iLow)-y(iLow-1))*(yLow - y(iLow)) ;
        else
            tLow = t(iLow) ;
        end
        yHigh = y(1) + RiseTimeLims(2)*(yfinal-y(1));
        iHigh = 1+ find((y(1:ns-1)-yHigh).*(y(2:ns)-yHigh) <= 0,1);
        if isempty(iHigh)
            // Has not reached the higher limit
            tHigh = %nan;
            settling_min = %nan ;
            settling_max = %nan ;
        else
            if (iHigh > 1) & (y(iHigh)~=y(iHigh-1))
                // Interpolate for more accuracy
                tHigh = t(iHigh) + (t(iHigh)-t(iHigh-1))/(y(iHigh)-y(iHigh-1))*(yHigh-y(iHigh));
            else
                tHigh = t(iHigh);
            end
            yRisen = y(iHigh:length(y)) ;
            settling_min = min(yRisen) ;
            settling_max = max(yRisen) ;
        end
        rise_time = tHigh - tLow ;
    else
        // System is not stable
        tHigh = %inf
        tLow = %inf
        rise_time = %nan ;
        settling_min = %nan ;
        settling_max = %nan ;
    end
endfunction

// Settling time with default 2% limits
// The user may change this by setting SettlingTimeThreshold to the desired value 
function[settlingTime] = settlingtime(y,t,yfinal,SettlingTimeThreshold)
    ns = length(t) ;
    if  abs(yfinal) < %inf
        err = abs(y-yfinal);
        rev_err = err(length(err):-1:1) ;
        tol = SettlingTimeThreshold* max(err);
        iSettle = length(err) - find(rev_err > tol, 1) + 1 ; 
        if isempty(iSettle)
            // Pure gain
            settlingTime = 0;
        elseif iSettle==ns
            //Has not settled
            settlingTime = %nan;
        elseif y(iSettle)~=y(iSettle+1)
            // Interpolate for more accuracy
            ySettle = yfinal + sign(y(iSettle)-yfinal) * tol;
            settlingTime = t(iSettle)+(t(iSettle)-t(iSettle+1))/(y(iSettle)-y(iSettle+1))*(ySettle-y(iSettle));
        else
            // Discrete time or pure gain
            settlingTime = t(iSettle+1);
        end 
    else
        // System is not stable
        settlingTime = %nan ;
    end
endfunction

// Main Step Info Function
function[RiseTime,SettlingTime,SettlingMin,SettlingMax,Overshoot,Undershoot,Peak,PeakTime,RiseTimeHigh,RiseTimeLow] = stepinfo(varargin)
    if(length(varargin) >= 3)
        for i = 1:length(varargin)-1
            if (typeof(varargin(i))=='string') & (varargin(i)=='RiseTimeLims') 
                if size(varargin(i+1))==[1 2]
                    RiseTimeLims = varargin(i+1) ;
                    if (size(RiseTimeLims)<>[1 2])|(RiseTimeLims(1)>RiseTimeLims(2))|(RiseTimeLims(2) > 1)|(RiseTimeLims(1) < 0)
                        error('The RiseTimeLims must be a real 1-by-2 vector with nondecreasing values between 0 and 1.') ;
                    end
                else
                    error('The RiseTimeLims must be a real 1-by-2 vector with nondecreasing values between 0 and 1.') ;
                end
            elseif (typeof(varargin(i))=='string') & (varargin(i)=='SettlingTimeThreshold') 
                SettlingTimeThreshold = varargin(i+1) ;
                if (SettlingTimeThreshold > 1) | (SettlingTimeThreshold < 0) 
                    error('SettlingTimeThreshold should be a scalar between 0 and 1') ;
                end
            else
                // Default values of RiseTimeLims and SettlingTimeThreshold
                RiseTimeLims = [0.1 0.9] ;
                SettlingTimeThreshold = 0.02 ;
            end
        end
    else
        // Default values of RiseTimeLims and SettlingTimeThreshold
        RiseTimeLims = [0.1 0.9] ;
        SettlingTimeThreshold = 0.02 ;
    end
    if typeof(varargin(1))== "constant" ;
        // Input is an array
        [ns ny nu] = size(varargin(1)) ;
        for k = 1:1:ny 
            for j = 1:1:nu
                y_temp = varargin(1) ;
                y = y_temp(1:ns,k,j)
                if (length(varargin) == 1)
                    // y_temp = varargin(1) ;
                    // y = y_temp(1:ns,k,j)
                    t = [1:1:length(y)] ;
                    yfinal =y(length(y)) ;
                elseif (length(varargin)==2)&(typeof(varargin(2)) == 'constant')
                    //y_temp = varargin(1) ;
                    //y = y_temp(1:ns,k,j) ;
                    t = varargin(2) ;
                    if (length(t) ~= length(y)) then
                        error("The input arguments Y and T must have compatible sizes.") ;
                    end
                    yfinal = y(length(y))
                elseif (length(varargin)==3)&(typeof(varargin(2)) == 'constant')&(typeof(varargin(3)) == 'constant') 
                    //y_temp = varargin(1) ;
                    //y = y_temp(1:ns,k,j) ;
                    t = varargin(2) ;
                    if length(t)~= length(y)
                        error('The input arguments Y and T must have compatible sizes.');
                    end
                    yfinal = varargin(3) ;
                    if size(yfinal) ~= [ny, nu]
                        error('The input arguments Y and YFINAL must have compatible sizes.');
                    end
                    yfinal = yfinal(k,j) ;
                elseif (length(varargin)==3)&(typeof(varargin(2)) == 'string')&(typeof(varargin(3)) == 'constant')
                    t = [1:1:length(y)] ;
                    yfinal =y(length(y)) ;
                elseif (length(varargin)>3)&(typeof(varargin(2))=='constant')&(typeof(varargin(3)) == 'string') 
                    t = varargin(2) ;
                    if (length(t) <> length(y)) then
                        error("The input arguments Y and T must have compatible sizes.") ;
                    end
                    yfinal = y(length(y))
                elseif (length(varargin)>3)&(typeof(varargin(2))=='constant')&(typeof(varargin(3)) == 'constant')
                    t = varargin(2) ;
                    if length(t)~= length(y)
                        error('The input arguments Y and T must have compatible sizes.');
                    end
                    yfinal = varargin(3) ;
                    if size(yfinal) ~= [ny, nu]
                        error('The input arguments Y and YFINAL must have compatible sizes.');
                    end
                    yfinal = yfinal(k,j) ;
                 else 
                     error('Wrong type of input arguments') ;
                end                
                [RiseTime(k,j), SettlingMin(k,j), SettlingMax(k,j),RiseTimeHigh(k,j),RiseTimeLow(k,j)]= risetime(y,t,yfinal,RiseTimeLims)
                SettlingTime(k,j) = settlingtime(y,t,yfinal,SettlingTimeThreshold)
                [Overshoot(k,j), Undershoot(k,j)] = overshoot_and_undershoot(y,t,yfinal) 
                [Peak(k,j), PeakTime(k,j)] = peak(y,t,yfinal) 
            end
        end
    elseif (typeof(varargin(1))== "rational" ) & (varargin(1).dt == 'c')
        // System is in rational form and continuous-time
        sys = varargin(1) ;
        [m n] = size(sys) ;
        for k = 1:1:m 
            for j = 1:1:n
                sys = varargin(1) ;
                sys = sys(k,j) ;
                if real(roots(sys.den))< 0
                    // system is stable
                    y_temp = csim('step',0:1:1010,sys) ;
                        yfinal = horner(sys,0.00001)
                        //Final Value theorem for calculation of yfinal
                        //In case of step input final value = limit s -> 0, C(z)
                    for i = [1:1:1000]  
                        temp = y_temp(i:i+9) ;
                        if (abs(temp- yfinal) < = 0.0001)
                            break
                            // Index of settling point
                        end
                    end 
                    // Simulate y only till the system settles
                    t = 0:(i+9)/10000:i+9 ;
                    y = csim('step',t,sys) ;
                else
                    // system is unstable
                    t = 0:(2.30*25)/((max(real(roots(sys.den)))+0.001)*10000):(2.30*25)/(max(real(roots(sys.den)))+0.001) ;
                    // t only till the value of y reaches around 10^25
                    y = csim('step',t,sys) ;
                    yfinal = sign(horner(sys,0.00001))*(%inf) ;
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit s -> 0, C(s)
                end
                [RiseTime(k,j), SettlingMin(k,j), SettlingMax(k,j),RiseTimeHigh(k,j),RiseTimeLow(k,j)]= risetime(y,t,yfinal,RiseTimeLims)
                SettlingTime(k,j) = settlingtime(y,t,yfinal,SettlingTimeThreshold)
                [Overshoot(k,j), Undershoot(k,j)] = overshoot_and_undershoot(y,t,yfinal) 
                [Peak(k,j), PeakTime(k,j)] = peak(y,t,yfinal) 
            end
        end
    elseif (typeof(varargin(1))=="state-space") & (varargin(1).dt == 'c')
        // System is in state-space form and continuous-time
        sys_ss = varargin(1) ;
        sys = ss2tf(sys_ss)
        // Converting to transfer function form
        [m n] = size(sys) ;
        for k = 1:1:m 
            for j = 1:1:n
                sys = ss2tf(sys_ss) ;
                sys = sys(k,j) ;
                if real(roots(sys.den))< 0
                    // system is stable
                    y_temp = csim('step',0:1:1010,sys) ;
                    yfinal = horner(sys,0.00001)
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit s -> 0, C(s)
                    for i = [1:1:1000]  
                        temp = y_temp(i:i+9) ;
                        if (abs(temp- yfinal) < = 0.0001)
                            break
                        end
                    end 
                    t = 0:(i+9)/10000:i+9 ;
                    y = csim('step',t,sys) ;
                else
                    // system is unstable
                    t = 0:(2.30*25)/((max(real(roots(sys.den)))+0.001)*10000):(2.30*25)/(max(real(roots(sys.den)))+0.001) ;
                    // t only till the value of y reaches around 10^25
                    y = csim('step',t,sys) ;
                    yfinal = sign(horner(sys,0.00001))*(%inf) ;
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit s -> 0, C(s)
                end
                [RiseTime(k,j), SettlingMin(k,j), SettlingMax(k,j),RiseTimeHigh(k,j),RiseTimeLow(k,j)]= risetime(y,t,yfinal,RiseTimeLims)
                SettlingTime(k,j) = settlingtime(y,t,yfinal,SettlingTimeThreshold)
                [Overshoot(k,j), Undershoot(k,j)] = overshoot_and_undershoot(y,t,yfinal) 
                [Peak(k,j), PeakTime(k,j)] = peak(y,t,yfinal) 
            end
        end
    elseif (typeof(varargin(1))=='rational') & (varargin(1).dt ~='c')
        // System is rational and discrete
        sys = varargin(1) ;
        [m n] = size(sys) ;
        for k = 1:1:m 
            for j = 1:1:n
                sys = varargin(1) ;
                sys = sys(k,j) ;
                //Calculation of sampling time
                if sys.dt == 'd'
                    sampling_time = 1 ;
                elseif sys.dt == []
                    error('Not a valid system') ;
                else
                    sampling_time = sys.dt ;
                end
                if abs(roots(sys.den))< 1
                    // system is stable
                    yfinal = horner(sys,1.0000)
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit z -> 1, C(z)
                    y = flts(ones(1,100000),sys) ;
                    t = [0:sampling_time:(length(y)-1)*sampling_time] ;
                else
                    // system is unstable
                    y = flts(ones(1,100000),sys) ;
                    t = [0:sampling_time:(length(y)-1)*sampling_time] ;
                    yfinal = sign(horner(sys,1.0000))*(%inf) ;
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit z -> 1, C(z)
                end
                [RiseTime(k,j), SettlingMin(k,j), SettlingMax(k,j),RiseTimeHigh(k,j),RiseTimeLow(k,j)]= risetime(y,t,yfinal,RiseTimeLims)
                SettlingTime(k,j) = settlingtime(y,t,yfinal,SettlingTimeThreshold)
                [Overshoot(k,j), Undershoot(k,j)] = overshoot_and_undershoot(y,t,yfinal) 
                [Peak(k,j), PeakTime(k,j)] = peak(y,t,yfinal) 
            end
        end
    elseif (typeof(varargin(1))=="state-space") & (varargin(1).dt ~= 'c')
        // System is in state-space form and discrete 
        sys_ss = varargin(1) ;
        sys = ss2tf(sys_ss)
        [m n] = size(sys) ;
        for k = 1:1:m 
            for j = 1:1:n
                sys = ss2tf(sys_ss) ;
                sys = sys(k,j) ;
                // Calculation of sampling time
                if sys.dt == 'd'
                    sampling_time = 1 ;
                elseif sys.dt == []
                    error('Not a valid system') ;
                else
                    sampling_time = sys.dt ;
                end
                if abs(roots(sys.den))< 1
                    // system is stable
                    yfinal = horner(sys,1.0000);
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit z -> 1, C(z)
                    y = flts(ones(1,100000),sys) ;
                    t = [0:sampling_time:(length(y)-1)*sampling_time] ;
                else
                    // system is unstable
                    y = flts(ones(1,100000),sys) ;
                    t = [0:sampling_time:(length(y)-1)*sampling_time] ;
                    yfinal = sign(horner(sys,1.0000))*(%inf) ;
                    //Final Value theorem for calculation of yfinal
                    //In case of step input final value = limit z -> 1, C(z)
                end
                [RiseTime(k,j), SettlingMin(k,j), SettlingMax(k,j),RiseTimeHigh(k,j),RiseTimeLow(k,j)]= risetime(y,t,yfinal,RiseTimeLims)
                SettlingTime(k,j) = settlingtime(y,t,yfinal,SettlingTimeThreshold)
                [Overshoot(k,j), Undershoot(k,j)] = overshoot_and_undershoot(y,t,yfinal) 
                [Peak(k,j), PeakTime(k,j)] = peak(y,t,yfinal) 
            end
        end
    else
        error('Wrong type of input argument') ;
    end
       disp(RiseTime,'RiseTime:') ;
       disp(SettlingTime,'SettlingTime:') ;
       disp(SettlingMin,'SettlingMin:') ;
       disp(SettlingMax,'SettlingMax:') ;
       disp(Overshoot,'Overshoot:') ;
       disp(Undershoot,'Undershoot:') ;
       disp(Peak,'Peak:') ;
       disp(PeakTime,'PeakTime:') ;
endfunction
    
