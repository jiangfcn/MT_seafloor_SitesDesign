function [linepar] = twopoints2line(x1,x2,y1,y2)
%line: ax-y+c = 0£»
atmp = (y1-y2)/(x1-x2);
btmp = y1-atmp*x1;
a = atmp;
b= -1;
c = btmp;
linepar = [a,b,c];
end
