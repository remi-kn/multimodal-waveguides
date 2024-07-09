function generate_geo_from_area_function(x, r, file_name, conical_segments)

    % Generates a csv file that can be used to import a geometry generated
    % from an area function in VocalTractLab3D
    
    % x                     x coordinate (cm)
    % r                     radius (cm)
    % file_name             name of the csv file to be generated
    %                       (Must end with '.csv')
    % conical_segments      if set to true, conical segments are generated
    %                       otherwise cylindrical segements are generated

    sep = ';';          % column separator for csv file
    n_x = length(x);
    
    if min(size(x)) == 1
        x = [x, zeros(size(x)), zeros(size(x))];
    end

    % generate circle contour
    nTheta = 120;
    theta = linspace(0, 2*pi, nTheta);
    cont = zeros(nTheta, 2);
    for ii = 1:nTheta
        cont(ii,1) = cos(theta(ii));
        cont(ii, 2) = sin(theta(ii));
    end

    no = [0. 1.];    % normal of all segments

    % compute scaling factors
    scaleIn = 1;
    scaleOut = ones(1,n_x);
    if conical_segments
        scaleOut(1:end-2) = 0.99*r(2:end-1)./r(1:end-2);
    end

    fid = fopen(file_name, 'w');

    for c = 1:n_x
        % write x coordinates
        fprintf(fid, '%f%s', x(c,1), sep); % x coordinate of center
        fprintf(fid, '%f%s', no(1), sep); % x coordinate of normal
        fprintf(fid, '%f%s', scaleIn, sep);    % entrance scaling factor
        % x coordinate of the contour
        for t = 1:nTheta
            fprintf(fid, '%f%s', r(c)*cont(t,1) + x(c,2), sep);  
        end
        fprintf(fid, '\n');

        % write y coordinates
        fprintf(fid, '%f%s', 0., sep); % y coordinate of center
        fprintf(fid, '%f%s', no(2), sep); % y coordinate of normal
        fprintf(fid, '%f%s', scaleOut(c), sep);    % exit scaling factor
        % y coordinate of the contour
        for t = 1:nTheta
            fprintf(fid, '%f%s', r(c)*cont(t,2) + x(c,3), sep);  
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end