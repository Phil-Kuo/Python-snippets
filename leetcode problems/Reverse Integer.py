class Solution:
    def reverse(self, x: int) -> int:
        ## n的符号
        sign = 1
        if x < 0:
            sign = -1
        reverse_list = []
        
        x = abs(x)
        while x:
            r = x % 10
            x //= 10
            reverse_list.append(r)
        
        l = len(reverse_list)
        s = 0
        for i in reverse_list:
            l -= 1
            s += i*10**l
        
        s = s*sign
        if (s >2**31-1 or s<-2**31):
            return 0
        
        return s

s = Solution()
print(s.reverse(-2**32))