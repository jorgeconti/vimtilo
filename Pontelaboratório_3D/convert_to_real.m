function [output] =  convert_to_real(input,type)

global maiorY menorY novoMinY novoMaxY maiorX menorX novoMinX novoMaxX

if(type==1)
    output=novoMinX+(( input- menorX)/( maiorX - menorX))*(novoMaxX-novoMinX);
elseif(type==2)
    output=novoMinY+(( input- menorY)/( maiorY - menorY))*(novoMaxY-novoMinY);
else
    output=0;
end