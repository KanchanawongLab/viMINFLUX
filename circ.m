function z = circ( x,y,D )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% D diameter of the circle
r=sqrt(x.^2+y.^2);
z=double(r<D/2);
z(r==D/2)=0;

end

