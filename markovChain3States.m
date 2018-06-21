function [channel] = markovChain3States(P,lengthChain)
channel = zeros(1,lengthChain); % 3-state Markov chain (output vector).
channel(1) = randi(3);
P1 = cumsum(P,2);
for i = 2:lengthChain
    event = randi(10000)/10000;
    if channel(1,i-1) == 1
        if event < P1(1,1) % No switch
            channel(1,i) = 1;
        elseif event > P1(1,2) % Switch to state 3
            channel(1,i) = 3;
        else
            channel(1,i) = 2; % Switch to state 2
        end
    elseif channel(1,i-1) == 2
        if event < P1(2,1) % Switch to state 1
            channel(1,i) = 1;
        elseif event > P1(2,2) % Switch to state 3
            channel(1,i) = 3;
        else % No switch
            channel(1,i) = 2;
        end
    elseif channel(1,i-1) == 3
        if event < P1(3,1) % Switch to state 1
            channel(1,i) = 1;
        elseif event > P1(3,2) % No switch
            channel(1,i) = 3;
        else % Switch to state 2
            channel(1,i) = 2;
        end
    end
end