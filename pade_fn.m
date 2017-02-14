function [a,varargout] = pade_fn(omega_n,f_n,N,varargin)
%%
%Louis-François Arsenault, September 2010
%
%[a] = pade_fn(omega_n,f_n,N)
%or
%[a,fw,omega] = pade_fn(omega_n,f_n,N,wmin,wmax,nbr_w,eta)
%
%This function implement the N points Padé coefficients for a function f(iwn).
%following the recursive algorithm with modification such that every iteration is
%kept normalized
%
%
%
%INPUT:
% omega_n  : the POSITIVE Matsubara frequencies entered as a line vector
% f_n      : the values of the function for positive omega_n, entered as 
%            a line vector
% N        : number of points of f(iwn) will be used. Must be < length(omega_n)
% wmin     : smallest real frequency to be calculated
% wmax     : biggest real frequency to be calculated
% nbr_w    : number or real freqencies to be calculated
% eta      : small imaginary part of the energy
%
%OUTPUT
% a        : Vector containing the a's
% fw       : Function in real frequencies
% omega    : real frequencies
%%

    %Test to see if the input vectors are line vectors
    [line,column] = size(omega_n);
    if line > 1
        error('omega_n must be a line vector')
    end
    [line,column] = size(f_n);
    if line > 1
        error('f_n must be a line vector')
    end

    %We only look at the first N points of the function
    omega_n = omega_n(1:N);
    f_n = f_n(1:N);

    %Contructing the g's
    %This matrix is artificial is the sense that the lower triangular part is
    %NAN always but that does not matter as we are only interested, at the end,
    %at the diagonal
    g(1,:) = f_n;
    g(2,:) = ( f_n(1) - f_n  )./(( i*omega_n - i*omega_n(1) ).*f_n);
    for k = 3:N
        g(k,:) = ( g(k-1,k-1) - g(k-1,:) )./( (i*omega_n - i*omega_n(k-1)).*g(k-1,:) );
    end
    
    %The a's is the diagonal of g
    a = transpose(diag(g));

    if nargin > 3
        wmin = varargin{1};
        wmax = varargin{2};
        nbr_w = varargin{3};
        eta = varargin{4};
        
        
        %The function is calculated in real frequencies 
        omega = wmin:(wmax-wmin)/(nbr_w-1):wmax;
        z = omega + i*eta;
        A = 1.0;
        B = 1.0 + a(2)*(z-i*omega_n(1)); 
        P = A./B;
        for c = 3:N
            A = 1.0 + a(c)*(z-i*omega_n(c-1))./A;
            B = 1.0 + a(c)*(z-i*omega_n(c-1))./B;
            P = P.*A./B;
        end
        fw = P*a(1);
            
       
        varargout{1} = fw;
        varargout{2} = omega;
    end
end