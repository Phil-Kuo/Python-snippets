import collections
class Solution:
    def threeSum(self, nums):
        nums = sorted(nums)
        l = len(nums)
        solution_list = set()
        for i in range(l-2):
            if i>0 and nums[i]==nums[i-1]:
                continue
            dict_nums = collections.defaultdict()
            dict_nums = {-nums[i]-n:j for j, n in enumerate(nums[i+1:])}
            for j, n in enumerate(nums[i+1:]):                               
                if n in dict_nums:
                    solution_list.update([(nums[i],n,-nums[i]-n)])
                
        return list(map(list, solution_list))

if __name__ =="__main__":
    S = Solution()
    print(S.threeSum([-1,0,1,2,-1,-4]))