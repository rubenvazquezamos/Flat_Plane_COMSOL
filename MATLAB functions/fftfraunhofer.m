function Ps = fftfraunhofer(theta,Rs,k,x)
% ======================================
%
%   PS = FFTFRAUNHOFER(THETA,RS,K,X)
%
%   FFTFRAUNHOFER calculates the far field PS using the Fraunhofer integral
%   of a surface with reflection RS(X) at wavenumber K and at angles THETA
%
%   Noé Jiménez, Salford, October 2016
%
%=======================================
na = length(x);
dx = x(2)-x(1);
Ps = zeros(size(theta));
for ia=1:na
    Ps = Ps+Rs(ia).*exp(1i*k*x(ia)*sin(theta))*dx; 
end

%Rs(ia), x(ia) and k are SCALARS. sin(theta) is a VECTOR

%dx is needed even though it doesn't affect the shape of the resultant
%polar pattern. If we don't multiply by dx then each value of Ps is quite
%large. dx infinitesimalises each value of Ps.

%notice Ps is not indexed by ia. This might seem weird at first, but the
%sin(theta) in the exp(term) ensures that the term "outputs" as a vector,
%giving the full contribution of a given point on the surface at all
%observation angles.

%Each step of the sum therefore is changing the value in every cell of Ps,
%by adding the contribution of each point on the surface at every angle.