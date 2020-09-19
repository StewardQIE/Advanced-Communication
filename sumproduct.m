function cer = sumproduct(iteration,ProbCrossover)
%% Use the belief propagation and sum-product to simulate the channle transmition
% Input: Max iteration and crossover probability
% Output: Simulate the communication system and CER
% Produced by Steward Qie 18/09/2020
%% Initialize the parameters
m = 3;
n = 2^m-1;
k = n-m;
M = 1000;
MatZero = [0 0 0];
Hmatrix = [0 1 1 1 0 1 0;1 0 0 1 1 1 0;1 1 0 1 0 0 1]; % H matrix

%% Generate message
message = randi([0,1],M,k);
codeword(:,1) = rem(message(:,1)+message(:,2)+message(:,3),2);
codeword(:,2) = rem(codeword(:,1)+message(:,1)+message(:,4),2);
codeword(:,3) = rem(codeword(:,2)+message(:,1)+message(:,3),2);
codeword(:,4:7) = message;
code1 = reshape(codeword,[],1);

%% Simulated the system and calculate CER
for i = 1:1:length(ProbCrossover)
    p_c = ProbCrossover(i);
    noise_code = bsc(codeword,p_c);
    sigma = p_c/2*7/4*7/4;
    for j = 1:1:M
        for l = 1:1:n
            pz(j,l) = 1/(1+exp(2*(2*noise_code(j,l)-1)/sigma)); 
            % prior probability if = 0
            po(j,l) = 1/(1+exp(-2*(2*noise_code(j,l)-1)/sigma)); 
            % prior probability if = 1
        end
    end
    for j = 1:1:M
        Imax = 1; % first iteration
        for a = 1:1:m
            for b = 1:1:n % message from variable node to ckpt node
                if Hmatrix(a,b) == 0
                    Qz(a,b) = 0;
                    Qo(a,b) = 0;
                else
                    Qz(a,b) = pz(j,b);
                    Qo(a,b) = po(j,b);
                end
            end
        end
        while Imax <= iteration
            Imax = Imax+1;
            for a = 1:1:m % message from ckpt node to variable node
                w = Hmatrix(a,:); % ckpt node 'a'
                index1 = find(w); % corresponding variable node
                for b = 1:1:n
                    if Hmatrix(a,b) == 0
                        Sz(a,b) = 0;
                        So(a,b) = 0;
                    else
                        index2 = index1;
                        x = find(index2==b);
                        index2(x) = []; % the variable node in the neighborhood
                        Sz(a,b) = 0;
                        So(a,b) = 0;
                        for x1 = 1:2^length(index2)
                            for x2 = length(index2):-1:1
                                if rem(x1-1,2^(-x2+length(index2)+1)) >= 2^(-x2+length(index2))
                                    v(x1,x2) = 1;
                                else
                                    v(x1,x2) = 0; % all possible values of the neighborhood variable node
                                end
                            end
                        end
                        for x1 = 1:2^length(index2)
                            ckpt = 0;
                            Pz = 1;
                            Po = 1; % initialization
                            for x2 = 1:1:length(index2)
                                ckpt = rem(v(x1,x2)+ckpt,2);
                            end
                            if ckpt == 0 % calculate the possibility of case one
                                for x2 = 1:1:length(index2)
                                    if v(x1,x2) == 0
                                        Pz = Pz*Qz(a,index2(x2));
                                    else
                                        Pz = Pz*Qo(a,index2(x2));
                                    end
                                end
                                Sz(a,b) = Sz(a,b)+Pz; % from ckpt node to variable node                               
                            else  % calculate the possibility of the other case
                                for x2 = 1:1:length(index2)
                                    if v(x1,x2) == 0
                                        Po = Po*Qz(a,index2(x2));
                                    else
                                        Po = Po*Qo(a,index2(x2));
                                    end
                                end
                                So(a,b) = So(a,b)+Po; % from ckpt node to variable node
                            end
                        end
                    end
                end
            end
            %calculate the updated message from variable node to ckpt node
            for a = 1:1:m
                for b = 1:1:n
                    if Hmatrix(a,b) == 0
                        Qz(a,b) = 0;
                        Qo(a,b) = 0;
                    else
                        w = Hmatrix(:,b);
                        index2 = find(w);
                        y = find(index2==a);
                        index2(y) = []; % ckpt node in the neighborhood
                        if index2
                            Qz(a,b) = 1;
                            Qo(a,b) = 1;
                            for l = 1:1:length(index2)
                                Qz(a,b) = Qz(a,b)*Sz(index2(l),b);
                                Qo(a,b) = Qo(a,b)*So(index2(l),b);
                            end
                            Qz(a,b) = Qz(a,b)*pz(j,b);
                            Qo(a,b) = Qo(a,b)*po(j,b);
                        else
                            Qz(a,b) = pz(j,b);
                            Qo(a,b) = po(j,b); %variable node to ckpt node
                        end
                    end
                end
            end
            %initialize Qz and Qo
            for a = 1:1:m
                for b = 1:1:n
                    if Qz(a,b) == 0
                        TempQz(a,b) = 0;
                        TempQo(a,b) = 0;
                    else
                        TempQz(a,b) = Qz(a,b)/(Qz(a,b)+Qo(a,b));
                        TempQo(a,b) = Qo(a,b)/(Qz(a,b)+Qo(a,b));
                    end
                end
            end
            Qz = TempQz;
            Qo = TempQo;
            % caculate the posterior probability of this iteration
            for b = 1:1:n
                PosPz(b) = 1;
                PosPo(b) = 1; %initialize posterior probability
                w = Hmatrix(:,b);
                index4 = find(w);
                for l = 1:1:length(index4)
                    PosPz(b) = PosPz(b)*Sz(index4(l),b);
                    PosPo(b) = PosPo(b)*So(index4(l),b);
                end
                PosPz(b) = PosPz(b)*pz(j,b);
                PosPo(b) = PosPo(b)*po(j,b);
            end
            % Posterior probability normalization
            for  b = 1:1:n
                TempPz(b) = PosPz(b)/(PosPz(b)+PosPo(b));
                TempPo(b) = PosPo(b)/(PosPz(b)+PosPo(b)); % reliability of each variable node
            end         
            for b = 1:1:n
                if TempPo(b) > 0.5 % decoded word determination
                    Decodeword(j,b) = 1;
                else
                    Decodeword(j,b) = 0;
                end
            end        
            z = Decodeword(j,:);
            sy = rem(z*(Hmatrix'),2); % Calculate syndrome
            if sy == MatZero % the decoded word is correct, stop iteration
                Imax = 150;
            end
        end
    end
    % remove ckpt bits, obtain the message and caculate the BER
    for a = 1:1:M
        for b = 4:1:n
            DecodeMes(a,b-3) = Decodeword(a,b);
        end
    end
    if p_c <= 0.5
        [number(i),cer(i)] = symerr(DecodeMes,message);
    else
        [number(i),cer(i)] = symerr(abs(1-DecodeMes),message);
    end
end

