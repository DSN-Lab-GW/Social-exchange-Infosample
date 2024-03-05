function NLL = mymodel_Optimal_ORA(pars, data)
% k       = pars(1);
% beta1   = pars(2);
% c       = pars(3);
% risk    = pars(4);
% alpha_0 = pars(5);
% beta_0  = pars(6);

k       = -1/(1+exp(-pars(1))); % alpha (transformed to be between zero and one);
beta1   = 10+100/(1+exp(-pars(3)));
c       = 5/(1+exp(-pars(3)));
risk    = (1/(1+exp(-pars(1)))-0.5)*2; 
% alpha_0 = 1+25/(1+exp(-pars(5)));
% beta_0  = 1+25/(1+exp(-pars(6)));
% risk = 0;
alpha_0 = 1;
beta_0  = 1; 

% NLL = 0;
m = 2;
T = 25;

        DeltaQ = computeDeltaQ_Optimal_ORA(T, m, c, risk, alpha_0, beta_0);
        
        trialidx   = find(data.red + data.green < T);
        
        thistime   = data.red(trialidx) + data.green(trialidx) + 1;
        thischoice = data.choice(trialidx);
        
        linearidx           = sub2ind(size(DeltaQ), data.green(trialidx) + 1, thistime);
        DeltaQ_vectorized   = DeltaQ(:);
        
        % Log likelihood
        prediction = 1./(1+exp(- thischoice .* (beta1 * (DeltaQ_vectorized(linearidx) - k)))); 
%         NLL = NLL - sum(log(prediction));
        NLL = sum(log(prediction+eps));

        
%     lowLimits =  [-1  10  0  -1  1   1];
%     highLimits = [ 0  inf 5   1  26 26]; 
