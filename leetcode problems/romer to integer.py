class Solution:
    def romanToInt(self, s: str) -> int:
        mapset = {'I':1,'V':5,'X':10,'L':50,'C':100,'D':500,'M':1000}
        if not s:
            return 0
        if len(s)==1:
            return mapset[s]
        
        l = len(s)
        nums = []
        for i in range(l-1):
            if mapset[s[i]] < mapset[s[i+1]]:
                nums.append(-mapset[s[i]])
            else:
                nums.append(mapset[s[i]])
        nums.append(mapset[s[-1]])
        return sum(nums)