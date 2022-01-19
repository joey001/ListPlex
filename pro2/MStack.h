namespace ListPlex
{
    template <typename T>
    struct Stack
    {
        T *members;
        int sz;
        int cap;

        Stack() : sz(0), members(nullptr) {}
        //Stack(T *_list, int _sz) : sz(_sz), members(list) {}
        //emty stack
        Stack(int _cap) : sz(0), cap(_cap)
        {
            members = new T[cap];            
        }
        Stack(const Stack &_stk) : sz(_stk.sz), cap(_stk.cap)
        {
            members = new T[cap];
            memcpy(members, _stk.members, sizeof(T) * sz);
        }
        ~Stack() { delete[] members; }

        inline void push(const T &v)
        {
            assert(sz < cap);
            members[sz++] = v;
        }
        inline void pop()
        {
            assert(sz > 0);
            sz--;
        }
        inline T top()
        {
            assert(sz > 0);
            return members[sz - 1];
        }

        template <typename Func>
        void for_each(Func &&f)
        {
            for (int i = 0; i < sz; i++)
            {
                f(members[i]);
            }
        }

        inline void clear(){
            sz = 0;
        }

        inline void resize(const int _cap){
            cap=_cap;
        }
    };
}