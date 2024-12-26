% Some figures. To be run with outputs from kuramoto.m 


%% Signals, parameter errors, synch errors, ...
for nr = 1:NR

    % Signals
    figure((nr-1)*4+1)
    clf
    mesh(tgrid,xgrid,real(u{nr}))
    view(0,90)
    title('master','Interpreter','latex')
    xlabel('$t$','Interpreter','latex')
    ylabel('$x$','Interpreter','latex')
    set(gca,'FontSize',18)
    if swa(1)
        figure((nr-1)*4+2)
        mesh(tgrid,xgrid,real(v{nr}))
        view(0,90)
        title('slave (fixed $\theta$)','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        ylabel('$x$','Interpreter','latex')
        set(gca,'FontSize',18)
    end %if
    if swa(2)
        figure((nr-1)*4+4)
        mesh(tgrid,xgrid,real(w{nr}))
        view(0,90)
        title('slave ($\hat\theta(t)$)','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        ylabel('$x$','Interpreter','latex')
        set(gca,'FontSize',18)
    end %if

    % normalised MSE vs. time
    figure((nr-1)*4+5)
    if swa(1)
        synchmse = mean( abs(u{nr}-v{nr}).^2 ) ./ mean(abs(u{nr}).^2);
        semilogy(tgrid,synchmse,'k-');
        xlabel('$t$','Interpreter','latex')
        ylabel('$\bar\mathcal{E}^2(t)$','Interpreter','latex')
        hold on
    end %if
    if swa(2)
        synchmse = mean( abs(u{nr}-w{nr}).^2 ) ./ mean(abs(u{nr}).^2);
        semilogy(tgrid,synchmse,'r-');
        xlabel('$t$','Interpreter','latex')
        ylabel('$\bar\mathcal{E}^2(t)$','Interpreter','latex')
        hold on
    end %if
    grid on;
    hold off;
    set(gca,'FontSize',18)

    % synchronisation error vs. time
    if swa(1)
        figure((nr-1)*4+6)
        real_err = real(u{nr}-v{nr});
        mesh(tgrid,xgrid,real_err)
        view(30,30)
        title('Synch. error','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        zlabel('$\mathcal{E}(t,x)$','Interpreter','latex')
        view(30,30)
        set(gca,'FontSize',18)
    end %if
    if swa(2)
        figure((nr-1)*4+8)
        real_err = real(u{nr}-w{nr});
        mesh(tgrid,xgrid,real_err)
        view(30,30)
        title('Synch. error','Interpreter','latex')
        xlabel('$t$','Interpreter','latex')
        zlabel('$\mathcal{E}(t,x)$','Interpreter','latex')
        view(30,30)
        set(gca,'FontSize',18)
    end %if

    % Parameter estimation (squared) error
    figure((nr-1)*4+9)
    true_par = [alpha beta gamma]';
    clf
    if swa(2)
        semilogy(tgrid,abs(true_par(1)-param{nr}(1,1:T)).^2,'k-');
        hold on
        semilogy(tgrid,abs(true_par(2)-param{nr}(2,1:T)).^2,'b--');
        semilogy(tgrid,abs(true_par(3)-param{nr}(3,1:T)).^2,'r-.');
        hold off
    end %if
    grid on
    hold off
    xlabel('$t$','Interpreter','latex');
    legend('$|\alpha-\hat\alpha|^2$','$|\beta-\hat\beta|^2$','$|\gamma-\hat\gamma|^2$','interpreter','latex');
    ylim([1e-20 10]);
    set(gca,'FontSize',18)

end %nr
