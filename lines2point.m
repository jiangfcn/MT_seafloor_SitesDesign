function pointxy = lines2point(a1,b1,c1,a2,b2,c2)
% pointxy = lines2point(a1,b1,c1,a2,b2,c2)
% a,b,and c represent a line equation a*x + b*y + c = 0
% pointxy include the coordinate of the intersection point of these two lines

k1 = -1*a1/b1; d1 = -1*c1/b1;
k2 = -1*a2/b2; d2 = -1*c2/b2;
px = (d2-d1)/(k1-k2);
py = k1*px + d1;
pointxy = [px,py];

end