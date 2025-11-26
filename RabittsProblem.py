
nandk = [35, 3]
def solve(nandk):
    n = nandk[0]
    k = nandk[1]
    if n < 3:
        return 1
    return solve((n - 1, k)) + k * solve((n - 2, k))

print(solve(nandk))