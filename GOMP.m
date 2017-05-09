function [x, Out] = GOMP(A,y,groups,varargin)
% This is a l21 minimization solver for group orthogonal matching pursuit:
%   minimiza ||x||_{2,1}
%   subject to y = Ax,
% -------------------------------------------------------------------------
% Author: Beier ZHU
%         Tsinghua University
% -------------------------------------------------------------------------
% 
% ==========================Required inputs================================
% A -- an m x n matrix
% y -- an m x 1 vector
% group -- an n-entry vector whose i-th entry is the group number of x_i,
%          the group index should start with 1, end with c, where c is the 
%          # of group
% ==========================Optional inputs================================
% 'maxIter' -- maximum number of iterations
% 'StopTolerance' -- stopping tolerance
% ===============================Outputs===================================
% x -- last iterate (hopefully an approximate solution)
% Out.iter -- # of iterations taken

    % Test for number of required parametres
    if (nargin-length(varargin)) ~= 3
        error('Wrong number of required parameters');
    end
    %----------------------------------------------------------------------
    % Set parameters to their defaults
    %----------------------------------------------------------------------
    opts.tol = 1e-4;
    opts.maxIter = 1e3;
    %----------------------------------------------------------------------
    % Set parameters to user specified values
    %----------------------------------------------------------------------
    if (rem(length(varargin),2)==1)
        error('Options should be given in pairs');
    else
        for i=1:2:(length(varargin)-1)
            switch upper(varargin{i})
                case 'STOPTOLERANCE'
                    opts.tol = varargin{i+1};
                case 'MAXITER'
                    opts.maxit = varargin{i+1};
                otherwise
                    error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
    [m_A, n_A] = size(A);
    n = length(groups);
    [m_y, n_y] = size(y);
    if n ~= n_A
        error(['group length must equal to the # of columns']);
    end
    if n_y ~= 1
        error(['y must be a column vector']);
    end
    if m_A ~= m_y
        error(['A and y must have same # of row'])
    end
    % ----------------------------------------------------------------------
    x = zeros(n,1);
    Omega = [];
    A_Omega = [];
    r = y;
    iter = 1;
    varepsilon = 1e-4;
    c = max(groups);
    for j = 1:c
        g{j} = find(groups == j);
    end
    while true
        l = A'*r;
        for j = 1:c 
            lg(j) = mean(abs(l(g{j})));
        end
        [~, max_index] = max(lg);
        group_index = g{max_index};
        Omega = [Omega group_index'];
        A_Omega = [A_Omega A(:,group_index')];
        x_k = A_Omega\y;
        r = y - A_Omega*x_k;
        iter = iter + 1;
        if (iter > opts.maxIter) || (norm(r) <= opts.tol) || (norm(A'*r, inf) <= varepsilon)
                break;
        end
    end
    x(Omega) = x_k;
    Out.iter = iter;
    Out.residual = norm(r);
    Out.energy = norm(A'*r, inf);
    % Out.tmp = Omega;
end