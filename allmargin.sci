//   Sys = ALLMARGIN(SYS) provides detailed information about the gain, phase,
//   and delay margins and the corresponding crossover frequencies of the
//   SISO open-loop model SYS.
//
//   The output S is a structure with the following fields:
//     * GMF: all -180 deg crossover frequencies in rad/TimeUnit
//     * GM: corresponding gain margins (g.m. = 1/G where G is the
//       gain at crossover)
//     * PMF: all 0 dB crossover frequencies (in rad/TimeUnit)
//     * PM: corresponding  phase margins (in degrees)
//     * DM, DMF: delay margins (in the units specified
//       in SYS.TimeUnit for continuous-time systems, and in multiples of
//       the sample time for discrete-time systems) and corresponding
//       critical frequencies
//     * Stable: 1 if stable, 0 if unstable, and NaN
//       if stability cannot be assessed (as in the case of most FRD systems)
//
//   Sys = ALLMARGIN(MAG,PHASE,W,TS) computes the stability margins from the
//   frequency response data W, MAG, PHASE and the sampling time TS. ALLMARGIN
//   expects gain values MAG in absolute units and phase values PHASE in 
//   degrees.
//   
//   EXAMPLE :
//   sys = syslin('c',(s +8)/(s^3+7*s^2+8*s+6))
//   a = allmargin(sys1)
//
//// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/
// http://in.mathworks.com/help/control/ ; 
// http://in.mathworks.com/help/control/ref/allmargin.html ;
// https://en.wikipedia.org/wiki/Control_systems;
//
//
// Author (s):
// Sanchit Gupta & Ashutosh Kumar Bhargava
//-----------------------------------------------------------------------------------------------------------------------//
function [output] = allmargin(varargin)
    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs == 2 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    elseif rhs >= 5 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    if rhs == 1 then
        sysData = varargin(1)
        select typeof(sysData)
        case "rational" then
        case "state-space" then
            sysData = ss2tf(sysData)
        else
            error(msprintf(gettext("Incompatible input argument.")))
        end
        sizeData = size(sysData)
        if (isequal(sizeData,[1 1])) == %F then
            error(msprintf(gettext("Input model must be SISO type.")))
        end
    elseif rhs == 3 | rhs == 4 then
        tempSize = size(varargin(1))
        if typeof(varargin(3)) == 'hypermat' | isequal(gsort(varargin(3),'c','i'),varargin(3)) == %F then 
            error(msprintf(gettext("frequency must be non-negative real valued vector and sorted in increasing order.")))
        elseif typeof(varargin(1)) == 'hypermat' then
            error(msprintf(gettext("mag must be non-negative real valued vector.")))
        elseif typeof(varargin(1)) <> 'constant' | typeof(varargin(2)) <> 'constant' | typeof(varargin(3)) <> 'constant' then
            error(msprintf(gettext("mag, phase, freq must be real valued vector with equal dimension")))
        elseif isequal(tempSize,size(varargin(2))) == %F | isequal(tempSize,size(varargin(3))) == %F then
            error(msprintf(gettext("mag, phase, freq must be real valued vector with equal dimension")))
        //elseif size(find(phasemag(varargin(1))<>0),"r") ~= 0 | size(find(phasemag(varargin(2))<>0),"r") ~= 0 | size(find(phasemag(varargin(3))<>0),"r") ~= 0  then
            //error(msprintf(gettext("mag, phase, freq must be real valued vector with equal dimension")))
        end
    end
//------------------------------------------------------------------------------------------------------------------//

    if rhs == 1 then
// Calculating phase margin
// Code is taken from p_margin
        eps=1.e-7;// threshold used for testing if complex numbers are real or pure imaginary
        h = sysData
        if h.dt=="c" then  //continuous time case
            w=poly(0,"w");
            niw=horner(h.num,%i*w);
            diw=horner(h.den,%i*w);
            // |n(iw)/d(iw)|=1 <-- (n(iw)*n(-iw))/(d(iw)*d(-iw))=1 <--  (n(iw)*n(-iw)) - (d(iw)*d(-iw))=0
            w=roots(real(niw*conj(niw)-diw*conj(diw)),"e");
            //select positive real roots
            ws=real(w(find((abs(imag(w))<eps)&(real(w)>0)))); //frequency points with unitary modulus
            if ws==[] then
                phm=[];
                fr=[];
                //return
            end
            f=horner(h,%i*ws);
        else  //discrete time case
            if h.dt=="d" then
                dt=1;
            else
                dt=h.dt;
            end
            // |h(e^(i*w*dt))|=1 <-- h(e^(i*w*dt))*h(e^(-i*w*dt))
            z=poly(0,varn(h.den));
            sm=simp_mode();
            simp_mode(%f);
            hh=h*horner(h,1/z)-1;
            simp_mode(sm);
            //find the numerator roots
            z=roots(hh.num,"e");
            z(abs(abs(z)-1)>eps)=[];// retain only roots with modulus equal to 1
            w=log(z)/(%i*dt);
            ws=real(w(abs(imag(w))<eps&real(w)>0)); //frequency points with unitary modulus
            if ws==[] then
                phm=%inf;
                fr=[];
                //return
            end
            f=horner(h,exp(%i*ws*dt));
        end
        phi=atand(imag(f),real(f));// phase of the frequency response (in [-180 180])
        //avoid near 0 negative phases that will give phm=180 instead of -180
        phi(phi>-1e-12&phi<0)=0;
        //compute the margins
        phm=pmodulo(phi,360)-180;
        //select the min value together with associated frequency in Hz
        frp=ws///(2*%pi);
    //---------------------------------------------------------------------------------------------------------------------------
    // calculatin phase margin
    // code is taken from g_margin
        epsr=1.e-7;//used for testing if complex numbers are real
        eps1=1.e-7;//used for testing if complex numbers have a modulus near 1
        epssing=1e-10; //used for testing if arguments are not singular points of h
        if h.dt=="c" then  //continuous time case
            // get s such as h(s)=h(-s) and s=iw
            s=%i*poly(0,"w");
            //compute h(s)-h(-s)=num/den
            num=imag(horner(h.num,s)*conj(horner(h.den,s)))
            den=real(horner(h.den,s)*conj(horner(h.den,s)))
            //necessary condition
            w=roots(num,"e");
            ws=real(w(abs(imag(w))<epsr&real(w)<=0)) //points where phase is -180°
    
            //remove nearly singular points
            ws(abs(horner(num,ws))>=epssing*abs(horner(den,ws)))=[]
            if ws==[] then gm=%inf,fr=[],return,end
            mingain=real(freq(h.num,h.den,%i*ws))
        else  //discrete time case
            if h.dt=="d" then dt=1,else dt=h.dt,end
            //get z such as h(z)=h(1/z) and z=e^(%i*w*dt)
            //form hh=h(z)-h(1/z)
            z=poly(0,varn(h.den));
            sm=simp_mode();simp_mode(%f);hh=h-horner(h,1/z);simp_mode(sm)
            //find the numerator roots
            z=roots(hh.num,"e");
            z(abs(abs(z)-1)>eps1)=[]// retain only roots with modulus equal to 1
    
            //remove nearly singular points
            z(abs(horner(hh.num,z))>=epssing*abs(horner(hh.den,z)))=[];
    
            w=log(z)/(%i*dt)
            ws=real(w(abs(imag(w))<epsr)) //points where phase is -180°
            if ws==[] then gm=%inf,fr=[],return,end
            mingain=real(horner(h,exp(%i*ws*dt)))
        end
    
        k=find(mingain<0)
        if k==[] then 
            gm=%inf;
            fr=[];
            //return
            end
        mingain=abs(mingain(k));
        ws=abs(ws(k))// select positive frequency
        gm=1/mingain//-20*log(mingain)/log(10) //tranform into Db
        frg=ws //transform in Hz
    //------------------------------------------------------------------------------------------------------------------//
        for ii = 1: length(phm)
            delayData(ii,1) = phm(ii,1)*2*%pi/(360*frp(ii,1))
            if delayData(ii,1) < 0 then
                delayData(ii,1) = delayData(ii,1)+2*%pi/frp(ii,1)
            end
        end
        if sysData.dt ~= 'c' then
            if sysData.dt == 'd' then
                sysData.dt = 1
            end
            delayData = abs(delayData)/sysData.dt
        end
    //-------------------------------------------------------------------------------------------------------------//
    // stability 
    stable = 0
    if sysData.dt == 'c'  
        // Continuous system
        if real(roots(sysData.den))< 0
            // system is stable
            stable = 1 ;
        end
    else
        //Discrete System
        if (real(roots(sysData.den))< 1)& (real(roots(sysData.den))> -1)
            // system is stable
            stable = 1 ;
        end
    end
    //--------------------------------------------------------------------------------------------------------------//
          output = struct('GM',gm,'GMF',frg','PM',phm','PMF',frp','DM',delayData','DMF',frp','stable',stable)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    elseif rhs == 3 | rhs == 4 then
        eps=1.e-2;
        magData = varargin(1)
        phaseData = varargin(2)
        freqData = varargin(3)
        nf = length(freqData);  
        freqData = matrix(freqData,[1 nf]);
        magData = matrix(magData,[1 nf]);
        phaseData = matrix(phaseData,[1 nf]);
        logmag = zeros(1,nf);
        isZero = (magData==0);
        logmag(:,isZero) = -%inf;
        logmag(:,~isZero) = log10(magData(~isZero));
        phaseData = unwrap((%pi/180)*phaseData);
        if nf>2 & freqData(1)==0
            freqData(1) = eps * freqData(2);  
        end
        logw  = log10(freqData);  
        twopi = 2*%pi;
        k = floor((phaseData+%pi)/twopi);
        lowcross = (2*k(:,1:nf-1)-1)*%pi;
        //lowcross = (2*k(:,1)-1)*%pi;
        ic = find(phaseData(:,2:nf)<lowcross | phaseData(:,2:nf)>=lowcross+twopi);
        //ic = ic(1:min(50,end)); 
        ic = ic(1:min(50,length(ic))); 
        Pc = lowcross(:,ic) + twopi*(phaseData(:,ic+1)>phaseData(:,ic));
        t = (Pc - phaseData(:,ic)) ./ (phaseData(:,ic+1) - phaseData(:,ic));
        frg = logw(:,ic) + t .* (logw(:,ic+1)-logw(:,ic));
        gm = logmag(:,ic) + t .* (logmag(:,ic+1)-logmag(:,ic));
        tol = %pi/6;
        if nf>=2,
            // Extrapolation toward freqData=0
            pcs = (2*round((phaseData(1)+%pi)/twopi)-1)*%pi; 
            if abs(phaseData(1)-pcs)<tol & abs(phaseData(2)-phaseData(1))~=0,  
                t = (pcs-phaseData(1)) / (phaseData(2)-phaseData(1));
                if t<0, 
                    frg = [frg , logw(1) + t * (logw(2)-logw(1))];
                    gm = [gm , logmag(1) + t * (logmag(2)-logmag(1))];
                end
            end
            // Extrapolation toward freqData=%inf
            pce = (2*round((phaseData(nf)+%pi)/twopi)-1)*%pi; 
            if abs(phaseData(nf)-pce)<tol & abs(phaseData(nf)-phaseData(nf-1))~=0, 
                t = (pcs-phaseData(nf-1)) / (phaseData(nf)-phaseData(nf-1));
                if t>0, 
                    frg = [frg , logw(nf-1) + t * (logw(nf)-logw(nf-1))];
                    gm = [gm , logmag(nf-1) + t * (logmag(nf)-logmag(nf-1))];
                end
            end
        end
        if isempty(gm)
            gm = zeros(1,0);   frg = zeros(1,0);
        else
            [frg,is] = gsort(frg);
            gm = 10.^(-gm(is));
            frg = 10.^frg;
        end
        // Phase margins  calculation(0dB gain crossings)
        ic = find(logmag(:,1:nf-1) .* logmag(:,2:nf) <= 0 & logmag(:,1:nf-1)~=logmag(:,2:nf));
        ic = ic(1:min(50,length(ic)));  
        t = -logmag(:,ic) ./ (logmag(:,ic+1) - logmag(:,ic));
        frp = logw(:,ic) + t .* (logw(:,ic+1)-logw(:,ic));
        phm = phaseData(:,ic) + t .* (phaseData(:,ic+1)-phaseData(:,ic));
        if isempty(phm)
            phm = zeros(1,0);   frp = zeros(1,0);
        else
            [frp,is] = gsort(frp);
            phm = pmodulo(phm(is),twopi)-%pi;
            frp = 10.^frp;
            phm = (180/%pi) * phm; 
        end
        // Delay Data
        for ii = 1: length(phm)
            delayData(ii,1) = phm(ii,1)*2*%pi/(360*frp(ii,1))
            if delayData(ii,1) < 0 then
                delayData(ii,1) = delayData(ii,1)+2*%pi/frp(ii,1)
            end
        end
        stable = %nan ;
        output = struct('GM',gm,'GMF',frg','PM',phm','PMF',frp','DM',delayData','DMF',frp','stable',stable)
    end        
endfunction






