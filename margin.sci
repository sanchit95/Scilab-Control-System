//   [Gm,Pm,Wcg,Wcp] = MARGIN(SYS) computes the gain margin Gm, the phase 
//   margin Pm, and the associated frequencies Wcg and Wcp, for the SISO 
//   open-loop model SYS (continuous or discrete). The gain margin Gm is 
//   defined as 1/G where G is the gain at the -180 phase crossing. The 
//   phase margin Pm is in degrees.
//
//  [Gm,Pm,Wcg,Wcp] = MARGIN(MAG,PHASE,W) derives the gain and phase margins 
//   from the Bode magnitude, phase, and frequency vectors MAG, PHASE, and W 
//   produced by BODE. MARGIN expects gain values MAG in absolute units and 
//   phase values PHASE in degrees. 
//   
//   MARGIN(SYS), by itself, plots the open-loop Bode plot with the gain 
//   and phase margins marked with a vertical line. 
//   EXAMPLE :
//   sys1 = syslin('c',(s +8)/(s^3+7*s^2+8*s+6))
//   [Gm,Pm,Wcg,Wcp] = allmargin(sys1)
//
//// References :
// http://www.scilab.org/resources/documentation ; 
// http://spoken-tutorial.org/
// http://in.mathworks.com/help/control/ ; 
// http://in.mathworks.com/help/control/ref/margin.html ;
// https://en.wikipedia.org/wiki/Control_systems;
//
//
// Author (s):
// Sanchit Gupta & Ashutosh Kumar Bhargava
//-----------------------------------------------------------------------------------------------------------------------//

function [gm,phm,frg,frp] = margin(varargin)
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
        if (lhs == 0)| (lhs == 1)
            show_margins(sysData) ;
        end    
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
            //Dm = utComputeDelayMargins(phm,frp,Ts,Td);
            phm = (180/%pi) * phm; 
        end
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //----Code Taken from show_margin() scilab function----//
        if (lhs == 0)|(lhs == 1)
            my_plt_margins(varargin(1),varargin(2),varargin(3),gm,phm,frg,frp) ;
        end
    end        
endfunction

// Code to show margins in bode plot
// Code taken from show_margins() scilab function and editted as required.
function my_plt_margins(magData,phaseData,freqData,gm,phm,frg,frp)
    fig=gcf();
    immediate_drawing=fig.immediate_drawing;
    fig.immediate_drawing="off";

    clf();
    f11 =matrix(freqData,[1 length(freqData)]) ; 
    p11 =matrix(phaseData,[1 length(phaseData)]) ;
    m11 = matrix(magData,[1 length(magData)])
    bode(f11,m11,p11) ;
    f=gcf();
    axg=f.children(2);
    axp=f.children(1);
    fmin=min(axg.x_ticks.locations);
    fmax=max(axg.x_ticks.locations);
    gmin=min(axg.y_ticks.locations);
    gmax=max(axg.y_ticks.locations);
    pmin=min(axp.y_ticks.locations);
    pmax=max(axp.y_ticks.locations);

    //[gm,fr]=g_margin(h)
    gm = max(gm)
    fr = max(frg)
    sca(axp);
    xpoly([fmin;fmax],[-180;-180])
    e=gce();e.foreground=color("red");e.line_style=4;
    if fr<>[] then
        xpoly([fr;fr],[pmin;pmax])
        e=gce();e.foreground=color("red");e.line_style=4;
        sca(axg);
        xpoly([fr;fr],[gmin;gmax])
        e=gce();e.foreground=color("red");e.line_style=4;
        xpoly([fr;fr],[-gm;0])
        e=gce();e.foreground=color("red");e.thickness=2;
    end
    //[phm,fr]=p_margin(h)
    phm = max(phm);
    fr = max(frp) ;
    sca(axg);
    xpoly([fmin;fmax],[0;0])
    e=gce();e.foreground=color("blue");e.line_style=4;
    if fr<>[] then
        xpoly([fr;fr],[gmin;gmax])
        e=gce();e.foreground=color("blue");e.line_style=4;
        sca(axp);
        xpoly([fr;fr],[pmin;pmax])
        e=gce();e.foreground=color("blue");e.line_style=4;
        xpoly([fr;fr],[-180;phm-180])
        e=gce();e.foreground=color("blue");e.thickness=2;
    end
    fig.immediate_drawing=immediate_drawing;
endfunction

//========================================================================================================================//







