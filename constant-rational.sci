// This is a sub function of series interconnection function.
//This contains the constant-rational combination
// The variable sys gives the resultant system ie. series(sys1,sys2,e,f)

function[sys] = series3(varargin)
    isSisoArray = %f ;
    for i =1:length(varargin) 
        if (typeof(varargin(i))== 'rational')| (typeof(varargin(i))== 'constant')|(typeof(varargin(i))== 'hypermat')
            sys1 = varargin(i)
            if (i+1<=length(varargin))& (typeof(varargin(i+1))=='boolean')&(varargin(i+1)==%T)
                isSisoArray = %T ;
            end
            break ;
        end
    end
    for i =1:length(varargin) 
        if ((typeof(varargin(i))== 'constant')|(typeof(varargin(i))== 'rational')|(typeof(varargin(i))== 'hypermat')) & (typeof(varargin(i))<>typeof(sys1))
            sys2 = varargin(i) 
            if (i+1<=length(varargin))&(typeof(varargin(i+1))=='boolean')&(varargin(i+1)==%T)
                isSisoArray = %T ;
            end
            break ;
        end
    end
    e = varargin(length(varargin)-1) ;  // vectored output
    f = varargin(length(varargin)) ;    // vectored input
    if (typeof(sys1)=='rational')
        //interchange sys1 and sys2
        temp = sys1;
        sys1 = sys2 ;
        sys2 = temp ; 
        //interchange values of e and f 
        temp = e;
        e = f ;
        f = temp ; 
    end
    if (typeof(sys1)=='constant')&(typeof(sys2)=='rational') & (ndims(sys2)<=2) &(isSisoArray == %f)
        // sys1 = constant type and sys2 = rational type 
        // Here sys2 is a pxq mimo system
        [m n] = size(sys1) ; 
        [p q] = size(sys2) ;  
        if ([m n]==[1 1]) & (e<=m)& (f<=q)
            // sys1 is 1x1 but sys2 is pxq 
            sys = sys1.*sys2(:,f) ;
        else 
            for i = 1:n
                for j= 1:p
                    sys(j,i) = sys1(e,i).*sys2(j,f) ; 
                end
            end
        end
    elseif (typeof(sys1)=='hypermat')&(typeof(sys2)=='rational')&(isSisoArray == %f)
        // Hypermatrix and Rational pxq mimo system
        ee=1;
        [p q] = size(sys2)
        for j=1:size(sys1,2)
            for k=1:size(sys1,3)
                temp(ee,1)=sys1(e,j,k)
                ee=ee+1
            end
        end
        for i =1:size(temp,1)
            for j = 1:size(sys2,1)
                sys(j,i) = temp(i,1).*sys2(j,f);
            end
        end
    elseif (typeof(sys2)=='rational')&(isSisoArray == %t) ;
        // sys2 is SISO array or a SISO hypermatrix of dimensions pxqxr 
        // Should produce error if f is not equal to 1
        p = size(sys2,1) ;
        q = size(sys2,2) ;
        r = size(sys2,3) ;
        ee = 1;
        temp = [] ;
            for j=1:size(sys1,2)
                for k=1:size(sys1,3)
                    //disp('loop');
                    temp(ee,1)=sys1(e,j,k);
                    ee=ee+1;
                end
            end
        temp_sys = [] ;
        for i =1:size(temp,1)
            temp_sys(i,:,:,:) = temp(i,1).*sys2;
        end
        sys = temp_sys ;
    end
endfunction
