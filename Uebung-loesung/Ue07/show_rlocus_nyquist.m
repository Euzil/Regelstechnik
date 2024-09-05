function []=show_rlocus_nyquist(sys)
    system = sys{1};
    name=sys{2};
    subplot(1,2,1)
    rlocus(system)
    subplot(1,2,2)
    nyquist(system)
    suptitle(name)
end
