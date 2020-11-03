function [x, e] = unitaleig_plot(num) % num +1 is the number of points

tic

x = zeros(1, num + 1); % grid over the MUB parameter
for k = 0 : num
  x(k + 1) = k * 2*pi / num;
end
e = zeros(1, num + 1); % maximal eigenvalue of the relevant operator

% Paulis

Id = eye(2);
X = [0 1; 1 0];
Y = [0 -i; i 0];
Z = [1 0; 0 -1];

% MUM operators
 
V = zeros(2, 2, 4);
V2 = zeros(2, 2, 4);
 
V(:, :, 1) = (X + Y + Z)/sqrt(3);
V(:, :, 2) = (X - Y - Z)/sqrt(3);
V(:, :, 3) = (- X + Y - Z)/sqrt(3);
V(:, :, 4) = (- X - Y + Z)/sqrt(3);
 
U = zeros(8);
P = zeros(8, 8, 4);
Q = zeros(8, 8, 4);
 
for j = 1 : 4,
 
    t = zeros(4);
    t(j, j) = 1;
 
    P(:, :, j) = kron( t, Id );
    U = U + kron( t, V(:, :, j) );
 
end
 
v = ones(4, 1)/2;
 
Q(:, :, 1) = kron( v * v', Id );
Q(:, :, 2) = U * Q(:, :, 1) * U';
 
X = 1/2*eye(8);
X(1, 4) = -1;
X(1, 5) = -1;
X(1, 8) = 1;
X(2, 3) = 1;
X(2, 6) = -1;
X(2, 7) = -1;
X(3, 6) = -1;
X(3, 7) = -1;
X(4, 5) = 1;
X(4, 8) = -1;
X(5, 8) = -1;
X(6, 7) = 1;
 
Q(:, :, 3) = (X + transpose(X))/4;
Q(:, :, 4) = eye(8) - Q(:, :, 1) - Q(:, :, 2) - Q(:, :, 3);

% MUBs

A = zeros(4, 4, 4);
B = zeros(4, 4, 4);

for a = 1 : 4
  A(a, a, a) = 1;
end

v2 = 1/2*[1; 1; -1; -1];
B(:, :, 1) = v * v';
B(:, :, 2) = v2 * v2';

p0 = [1, 2, 3, 4]; % for the post-processing

% Loop over the MUB parameter
for k = 1 : num + 1
  % construct the x-dependent MUB vectors and projectors
  v3 = 1/2 * [1; -1; 1i * exp( x(k) * 1i ); -1i * exp( x(k) * 1i )];
  v4 = 1/2 * [1; -1; -1i * exp( x(k) * 1i ); 1i * exp( x(k) * 1i )];
  B(:, :, 3) = v3 * v3';
  B(:, :, 4) = v4 * v4';
  % loop over post-processings
  for m = 1 : 24
    p = perms(p0)(m, :); % set the permutation
    S = 0; % to store the relevant operator

    for a = 1 : 4
      S = S + kron(transpose( P(:, :, a) ), A(:, :, a) ) +... 
      kron(transpose( Q(:, :, p(a)) ), B(:, :, a) );
    end
    e_tmp = max(real(eig(S))); % current value of the maximal eigenvalue

    if e_tmp >= e(k)
      e(k) = e_tmp; % maximal eigenvalue over permutations
    end
  end
end

toc

end