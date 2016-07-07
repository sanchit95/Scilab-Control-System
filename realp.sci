// Function realp()
// [variable, details] = realp(NAME,VALUE)
// creates a tunable real-valued parameter 'x' with name and initial value specified by the string NAME and the numeric array 
// VALUE.
//Using ordinary arithmetic operators (+, -, *, /, \, ^), you can combine real parameter objects into rational expressions and 
//use these expressions to create parametric models (both static and dynamic). You can then use such models to perform 
//parameter studies, or tune control systems.
//
// NOTE :-
// Load PIMS loader in scilab to run this code, as this code depends on python sympy library and python varsion 2.7
// Link to PIMS Library for downloading and Installing
// https://atoms.scilab.org/toolboxes/PIMS
//
// EXAMPLE
// [x x_details] = realp('x',5)
// [y y_details] = realp('y',3) 
//  eqn = x^y + 5*x*y + 89/x ;
// 
// Author (s):
// Sanchit Gupta & Ashutosh Kumar Bhargava
//----------------------------------------------------------------------------------------------------------------------//
function [variable,details] = realp(varargin)
    pyImport sympy
    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=3 | rhs == 1 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    varbName = varargin(1)
    varbData = varargin(2)
    if typeof(varargin(1)) <> 'string' then
        error((gettext("Wrong type for argument: String expected.")))//msprintf
    elseif typeof(varargin(2)) <> 'constant' then
        error((gettext("Wrong type for argument: Real matrix (2D) expected.")))//msprintf
    end
    variable = sympy.symbols(varargin(1))
    varbPhase = phasemag(varbData)
    findImag = find(varbPhase <> 0 & varbPhase <> 180)
    if size(findImag,"r") ~= 0 then
        error((gettext("Wrong type for argument: Real matrix (2D) expected.")))
    end
    sizeData = size(varbData)
    numbOfElements = 1 
    for ii = 1 : 2
        numbOfElements = numbOfElements*sizeData(ii)
    end
    for ii = 1: numbOfElements
        minData(ii,1) = -%inf
        freeData(ii,1) = 1
    end
    maxData = -minData
    minData = hypermat([sizeData],minData)
    maxData = hypermat([sizeData],maxData)
    freeData = hypermat([sizeData],freeData)
    if (size(minData)<>size(maxData))
        error('Cannot change the size of realp data') ;
    end
    details("Name") = varbName
    details("Value") = varbData
    details("Minimum") = minData
    details("Maximum") = maxData
    details("Free") = freeData
endfunction
