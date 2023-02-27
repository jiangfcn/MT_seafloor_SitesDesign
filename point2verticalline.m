function [linepar] = point2verticalline(x,y,a,b)
%[linepar] = point2verticalline(x,y,a,b)
% vertical line: linepar(1)*x + linepar(2)*y + linepar(3) = 0£»
%passing a known point to draw a vertical line of a known line (ax+by+c=0);
% and solving the equation of the vertical line(y = a2x+c2)

k = b/a;%slop of a vertical line.
a2 = k;
b2 = -1;
c2 = -1*k*x+y;
linepar = [a2,b2,c2]; 
end