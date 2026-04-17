function F = hankel_transform(Fr, r, dist)
    % More robust Hankel transform implementation
    F = zeros(size(dist));

    for i = 1:length(dist)
        if dist(i) == 0
            % Special case for k=0
            integrand = Fr ;
        else
            integrand = Fr .* besselj(0, dist(i) * r);
        end
        
        % Use trapezoidal rule for integration
        F(i) = 1/(2*pi) * trapz(r, integrand);
        if i==length(dist)/2
            fprintf('half way there... ');
        end
    end
end
