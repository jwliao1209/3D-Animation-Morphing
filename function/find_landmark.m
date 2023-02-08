function ID = find_landmark(V, P)
    [~,ID] = min(sum((V-P).^2, 2));
end