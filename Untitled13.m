% PROGRAM NAME: ps4huggett.m
clear, clc

%Problem 1
% define FE: v(s,a)=max U(y(s)+a-qa')+beta*Es'|s[v(s',a')]
%policy function choose a' from constraint correspondence
% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix



% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 10;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector
a_mat = repmat(a', [1 num_a]);
% INITIAL GUESS FOR q, q is constant , try diffrent q if nessasery
q_min = 0.98;
q_max = 1;
%q_guess = (q_min + q_max) / 2;
q_guess =0.9985;
% ITERATE OVER ASSET PRICES
%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
Y_s=zeros(num_a,num_a,2);
for i=1:num_a
    Y_s(i,:,1)=1;
    Y_s(i,:,2)=b;
end 

dis1 = 1 ;

cons(:,:,1) =  Y_s(:,:,1)+  a_mat - q_guess*a_mat'; 
% from the constraint, in the A high case, generate all c given different
% combo of k and k'
cons(:,:,2) = Y_s(:,:,2) +  a_mat - q_guess*a_mat';

ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
% given k and k', we can get the No. value of the value funtion
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
% v_guess = zeros(2, num_k);
v_guess = zeros(1, num_a,2);
% set initial value of value function to 0
while dis > tol
    % an alternative, more direct way to compute the value array:
    % value_mat _alt = ret + beta * ...
       % repmat(permute((T_mat * v_guess), [3 2 1]), [num_k 1 1]);
    
    % compute the utility value for all possible combinations of k and k':
   % finish from here!
    value_mat(:,:,1) = ret(:,:,1) + beta * ( ...
        PI(1,1) * repmat(v_guess(1,:,1), [num_a 1]) + ...
        PI(1,2) * repmat(v_guess(1,:,2), [num_a 1])); 
   
    value_mat(:,:,2) = ret(:,:,2) + beta * ( ...
       PI(2,1) * repmat(v_guess(1,:,1), [num_a 1]) + ...
        PI(2,2) * repmat(v_guess(1,:,2), [num_a 1]));
    
    % find the optimal k' for every k:
    [vfn, pol_indx] = max(value_mat, [], 2);
    % find max value for each row
    % vfn = vfn';
    vfn=permute(vfn,[2 1 3]);
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

    
    % KEEP DECSISION RULE
    pol_fn = a(pol_indx);   
    
    % SET UP INITITAL DISTRIBUTION
 % MASS(1,:)=mass(1:num_a,1);
   % MASS(2,:)=mass(num_a+1:2*num_a,1);
    Mu_guess=zeros(num_a,2);
    Mu_guess(:,:)=1/(2*num_a);
  
    % ITERATE OVER DISTRIBUTIONS
   while dis1>0.0001;
    [emp_ind,a_ind,mass] = find(Mu_guess > 0); % find non-zero indices
  
    %emp_ind=employment index, a_ind=assest index, 
    MuNew = zeros(size(Mu_guess));
   % give same dim for MuNew
  %for ii = 1:length(emp_ind)
   %     apr_ind = pol_indx(a_ind(ii),1,emp_ind(ii));
    %    MuNew(:, apr_ind) = MuNew(:, apr_ind) + (PI(emp_ind(ii), :) * mass)';
   %end 
   for i=1:num_a;
       k=find(pol_indx==i);
       k1=k(k<11);
       k2=k(k>10);
   MuNew(i,1)=sum(PI(1,1)*Mu_guess(k1,1))+sum(PI(2,1)*Mu_guess(k2-num_a,2));
   MuNew(i,2)=sum(PI(1,2)*Mu_guess(k1,1))+sum(PI(2,2)*Mu_guess(k2-num_a,2));
   
   end 
   dis1 = max(abs(MuNew - Mu_guess));
   Mu_guess=MuNew;

   end 
   %above is sln to problem 2
  for i=1:num_a
   Asset(i)=pol_fn(i,1,1)*MuNew(i,1)+pol_fn(i,1,2)*MuNew(i,2);
  end 
Total=sum(Asset);
% adjust q to let Total = 0
   
 
        
